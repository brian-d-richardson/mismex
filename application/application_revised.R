# Load libraries
rm(list = ls())
library(geex)
library(ggplot2)
library(mismex)
library(MASS)
library(dplyr)

# Set seed for full analysis
set.seed(1234)

#########################################################################
# Data Processing
#########################################################################

# Build application data set
load("HVTN505/data/dat.505.rda")

assays <- subset(var.505, assay %in% c("fcrR2a", "fcrR3a", "phago"))

primarydat <- read.csv("primary505_for_sharing_upd.csv")
primarydat$ptid <- primarydat$pub_id

fulldat <- merge(dat.505, primarydat, by = "ptid", all = T)
fulldat$marker1 <- log(fulldat$ADCP1)
fulldat$marker2 <- fulldat$R2aConSgp140CFI

# Treat behavior risk as binary variable
fulldat$bhv_bin.y <- 1*(fulldat$bhvrisk.y > 0)

# Simple imputation of missing confounders, marker2 0 value
bmimod <- lm(BMI.y ~ HIVall.y + age.y + bhvrisk.y + race.y, data = fulldat)
bmi_imp_dat <- data.frame(HIVall.y = fulldat$HIVall.y[is.na(fulldat$BMI.y)],
                          age.y = fulldat$age.y[is.na(fulldat$BMI.y)],
                          bhvrisk.y = fulldat$bhvrisk.y[is.na(fulldat$BMI.y)],
                          race.y = fulldat$race.y[is.na(fulldat$BMI.y)])
fulldat$BMI.y[is.na(fulldat$BMI.y)] <- predict(bmimod, bmi_imp_dat)

bhvmod <- glm(bhv_bin.y ~ HIVall.y + age.y + BMI.y + race.y, data = fulldat,
              family = binomial)
bhv_imp_dat <-
  data.frame(HIVall.y = fulldat$HIVall.y[is.na(fulldat$bhv_bin.y)],
             age.y = fulldat$age.y[is.na(fulldat$bhv_bin.y)],
             BMI.y = fulldat$BMI.y[is.na(fulldat$bhv_bin.y)],
             race.y = fulldat$race.y[is.na(fulldat$bhv_bin.y)])
fulldat$bhv_bin.y[is.na(fulldat$bhv_bin.y)] <-
  1*(predict(bhvmod, bhv_imp_dat, type = "response") > 0.5)

imputemax <- min(fulldat$marker2[fulldat$marker2 > 0], na.rm = T)
fulldat$marker2[fulldat$marker2 == 0 & !is.na(fulldat$marker2) &
                  fulldat$trt.y == 1] <- runif(1, 0, imputemax)

fulldat$sampled <- !is.na(fulldat$marker1)
sampmod <- glm(sampled ~ HIVwk28preunbl.y,
               data = fulldat[fulldat$trt.y == 1, ])

# Create analysis dataset
# Restrict to treatment, use both markers
fulldat$CD4PFS <- fulldat$CD4_ANYVRCENV_PolyfunctionalityScore_score_bin
fulldat$CD8PFS <- fulldat$CD8_ANYVRCENV_PolyfunctionalityScore_score_bin
analysisdat <- fulldat %>%
  filter(trt.y == 1) %>%
  select(HIVall.y, HIVwk28preunbl.y, marker1, marker2,
         race.y, BMI.y, age.y, bhv_bin.y, ADCP1_bin, CD4PFS, CD8PFS, ptid)

# Probability selected by treatment
# Weights based on HIV at week 28, not end of study
analysisdat$sampweights <-
  (!is.na(analysisdat$marker1))*(analysisdat$HIVwk28preunbl.y == 1) /
  (25 / sum(fulldat$HIVwk28preunbl.y[fulldat$trt.y == 1])) +
  (!is.na(analysisdat$marker1))*(analysisdat$HIVwk28preunbl.y == 0) /
  (125 / sum(1 - fulldat$HIVwk28preunbl.y[fulldat$trt.y == 1]))


#########################################################################
# Analysis for the first marker: ADCP
#########################################################################

# Remove individuals without the marker measured
analysisdat1 <- analysisdat %>%
  filter(!is.na(marker1))

# Define inverse identity link and its derivative
inv.ident <- function(x) { x }
d.inv.ident <- function(x) { rep(1, length(x)) }

# Define inverse log link and its derivative
inv.log <- function(x) { exp(x) }
d.inv.log <- function(x) { exp(x) }

# Subset data and rename variables
analysisdat1 <- analysisdat1[, c(2, 3, 5:8, 10:13)]
colnames(analysisdat1) <- c("Y", "A", "L1", "L2", "L3",
                            "L4", "L5", "L6", "ID", "cc.wts")

# Arguments for DR estimation
dr.args <- list(formula = "~ A + L1 + L2 + L3 + L4 + L5 + L6",
                ps.formula = "~L1 + L2 + L3 + L4 + L5 + L6",
                inv.link = inv.log,
                d.inv.link = d.inv.ident)

# Sequence of exposure values
m1vals <- seq(0.5, 3, 0.1)

# Measurement error values considered
me_list <- c(0, 1/8, 3/16, 1/4)

# Initialize DR estimates and confidence intervals
dr_ests1 <- array(NA, dim = c(4, length(m1vals)))
dr_lower1 <- array(NA, dim = c(4, length(m1vals)))
dr_upper1 <- array(NA, dim = c(4, length(m1vals)))

# Fit DR model with no measurement error
dr.res0 <- fit.dr(data = analysisdat1, args = dr.args, a = m1vals)

# Grab results corresponding to counterfactual curve
EYa.ind <- grepl("EYa", names(dr.res0$est))
EYa <- dr.res0$est[EYa.ind]
cov.EYa <- dr.res0$bc.var[EYa.ind, EYa.ind]

# Delta method to log
grad.log <- diag(1 / EYa)
se.log.EYa <- sqrt(diag(grad.log %*% cov.EYa %*% t(grad.log)))

# Add results with no measurement error
dr_ests1[1, ] <- EYa
dr_lower1[1, ] <- exp(log(EYa) - qnorm(0.975) * se.log.EYa)
dr_upper1[1, ] <- exp(log(EYa) + qnorm(0.975) * se.log.EYa)

# Results for increasing measurement error
for (me in 2:4) {

  sigma_me1 <- var(analysisdat1$A)*me_list[me]

  # DR estimator
  dr.res <- fit.dr.mccs(data = analysisdat1, a = m1vals,
                        cov.e = sigma_me1, B = 100, mc.seed = 1234,
                        return.var = TRUE, args = dr.args,
                        start = dr_ests1[me - 1, 1:8])

  # Grab results corresponding to counterfactual curve
  EYa.ind <- grepl("EYa", names(dr.res$est))
  EYa <- dr.res$est[EYa.ind]
  cov.EYa <- dr.res$bc.var[EYa.ind, EYa.ind]

  # Delta method to log
  grad.log <- diag(1 / EYa)
  se.log.EYa <- sqrt(diag(grad.log %*% cov.EYa %*% t(grad.log)))
  dr_ests1[me, ] <- EYa
  dr_lower1[me, ] <- exp(log(EYa) - qnorm(0.975) * se.log.EYa)
  dr_upper1[me, ] <- exp(log(EYa) + qnorm(0.975) * se.log.EYa)

}


#########################################################################
# Analysis for the second marker: RII
#########################################################################

# Remove individuals without the marker measured
analysisdat2 <- analysisdat %>%
  filter(!is.na(marker2))

# Subset data and rename variables
analysisdat2 <- analysisdat2[, c(2, 4:8, 10:13)]
colnames(analysisdat2) <- c("Y", "A", "L1", "L2", "L3",
                            "L4", "L5", "L6", "ID", "cc.wts")

# Arguments for DR estimation
dr.args2 <- list(formula = "~ A + L1 + L2 + L3 + L4 + L5 + L6",
                 ps.formula = "~L1 + L2 + L3 + L4 + L5 + L6",
                 inv.link = inv.log,
                 d.inv.link = d.inv.ident)

# Sequence of exposure values
m2vals <- seq(7, 10, 0.1)

# Measurement error values considered
me_list2 <- c(0, 1/8, 3/16, 1/4)

# Initialize DR estimates and confidence intervals
dr_ests2 <- array(NA, dim = c(4, length(m2vals)))
dr_lower2 <- array(NA, dim = c(4, length(m2vals)))
dr_upper2 <- array(NA, dim = c(4, length(m2vals)))

# Fit DR model with no measurement error
dr.res2.0 <- fit.dr(data = analysisdat2, args = dr.args2, a = m2vals)

# Grab results corresponding to counterfactual curve
EYa2.ind <- grepl("EYa", names(dr.res2.0$est))
EYa2 <- dr.res2.0$est[EYa2.ind]
cov.EYa2 <- dr.res2.0$bc.var[EYa2.ind, EYa2.ind]

# Delta method to log
grad.log2 <- diag(1 / EYa2)
se.log.EYa2 <- sqrt(diag(grad.log2 %*% cov.EYa2 %*% t(grad.log2)))

# Add results with no measurement error
dr_ests2[1, ] <- EYa2
dr_lower2[1, ] <- exp(log(EYa2) - qnorm(0.975) * se.log.EYa2)
dr_upper2[1, ] <- exp(log(EYa2) + qnorm(0.975) * se.log.EYa2)

# Results for increasing measurement error
for (me in 2:4) {

  sigma_me2 <- var(analysisdat2$A)*me_list2[me]

  # DR estimator
  dr.res2 <- fit.dr.mccs(data = analysisdat2, a = m2vals,
                         cov.e = sigma_me2, B = 100, mc.seed = 12345,
                         return.var = TRUE, args = dr.args2,
                         start = dr_ests2[me - 1, 1:8])

  # Grab results corresponding to counterfactual curve
  EYa2.ind <- grepl("EYa", names(dr.res2$est))
  EYa2 <- dr.res2$est[EYa2.ind]
  cov.EYa2 <- dr.res2$bc.var[EYa2.ind, EYa2.ind]

  # Delta method to log
  grad.log2 <- diag(1 / EYa2)
  se.log.EYa2 <- sqrt(diag(grad.log2 %*% cov.EYa2 %*% t(grad.log2)))
  dr_ests2[me, ] <- EYa2
  dr_lower2[me, ] <- exp(log(EYa2) - qnorm(0.975) * se.log.EYa2)
  dr_upper2[me, ] <- exp(log(EYa2) + qnorm(0.975) * se.log.EYa2)

}


#########################################################################
# Plot results
#########################################################################
latdat <- data.frame(vals = c(rep(m1vals, 4),
                              rep(m2vals, 4)),
                     Risk = c(dr_ests1[1, ], dr_ests1[2, ],
                              dr_ests1[3, ], dr_ests1[4, ],
                              dr_ests2[1, ], dr_ests2[2, ],
                              dr_ests2[3, ], dr_ests2[4, ]),
                     Risk_low = c(dr_lower1[1, ], dr_lower1[2, ],
                                  dr_lower1[3, ], dr_lower1[4, ],
                                  dr_lower2[1, ], dr_lower2[2, ],
                                  dr_lower2[3, ], dr_lower2[4, ]),
                     Risk_upp = c(dr_upper1[1, ], dr_upper1[2, ],
                                  dr_upper1[3, ], dr_upper1[4, ],
                                  dr_upper2[1, ], dr_upper2[2, ],
                                  dr_upper2[3, ], dr_upper2[4, ]),
                     ME = c(rep(round(me_list, 3), each = length(m1vals)),
                            rep(round(me_list2, 3), each = length(m2vals))),
                     Exposure = c(rep("ADCP", 4*length(m1vals)),
                                  rep("RII", 4*length(m2vals))))

# New facet label names for Exposure variable
ggplot(latdat, aes(x = vals, y = Risk, ymin = Risk_low, ymax = Risk_upp)) +
  geom_line() +
  facet_grid(ME ~ Exposure, scales = "free",
             labeller = label_bquote(sigma[me]^2 == .(ME) * sigma^2)) +
  #labeller = labeller(ME = me.labs, Exposure = exp.labs)) +
  geom_ribbon(alpha = 0.3) +
  xlab("Exposure values") + ylab("HIV risk") +
  ylim(c(0, 0.4)) +
  theme_bw()


#########################################################################
# Diagnostics
#########################################################################

# First will fit the component models outside of the new estimation functions
# ADCP propensity model
prop_mod1 <- lm(A ~ L1 + L2 + L3 + L4 + L5 + L6, data = analysisdat1)

# RII propensity model
prop_mod2 <- lm(A ~ L1 + L2 + L3 + L4 + L5 + L6, data = analysisdat2)

# Outcome model conditional on ADCP
out_mod1 <- fit.glm(analysisdat1, dr.args)
out_preds1 <- exp(out_mod1$est %*% t(cbind(1, analysisdat1[, 2:8])))

# Outcome model conditional on RII
out_mod2 <- fit.glm(analysisdat2, dr.args2)
out_preds2 <- exp(out_mod2$est %*% t(cbind(1, analysisdat2[, 2:8])))

# Propensity model diagnostics
par(mfrow = c(2, 2))
plot(prop_mod1)
plot(prop_mod2)

# Outcome model diagnostics (manual Chi-square GoF)
# ADCP outcome model
pi1_sum1 <- sum(sort(out_preds1)[1:38]*
                analysisdat1$cc.wts[order(out_preds1)[1:38]])
pi1_sum2 <- sum(sort(out_preds1)[39:75]*
                analysisdat1$cc.wts[order(out_preds1)[39:75]])
pi1_sum3 <- sum(sort(out_preds1)[76:112]*
                analysisdat1$cc.wts[order(out_preds1)[76:112]])
pi1_sum4 <- sum(sort(out_preds1)[113:150]*
                analysisdat1$cc.wts[order(out_preds1)[113:150]])

y1_sum1 <- sum(analysisdat1$Y[order(out_preds1)[1:38]]*
               analysisdat1$cc.wts[order(out_preds1)[1:38]])
y1_sum2 <- sum(analysisdat1$Y[order(out_preds1)[39:75]]*
               analysisdat1$cc.wts[order(out_preds1)[39:75]])
y1_sum3 <- sum(analysisdat1$Y[order(out_preds1)[76:112]]*
               analysisdat1$cc.wts[order(out_preds1)[76:112]])
y1_sum4 <- sum(analysisdat1$Y[order(out_preds1)[113:148]]*
               analysisdat1$cc.wts[order(out_preds1)[113:148]])

ts <- (y1_sum1 - pi1_sum1)^2 / (pi1_sum1*(1 - pi1_sum1/38)) +
  (y1_sum2 - pi1_sum2)^2 / (pi1_sum2*(1 - pi1_sum2/37)) +
  (y1_sum3 - pi1_sum3)^2 / (pi1_sum3*(1 - pi1_sum3/37)) +
  (y1_sum4 - pi1_sum4)^2 / (pi1_sum4*(1 - pi1_sum4/38))

pchisq(ts, 2, lower.tail = FALSE)

# RII outcome model
pi2_sum1 <- sum(sort(out_preds2)[1:38]*
                  analysisdat2$cc.wts[order(out_preds2)[1:38]])
pi2_sum2 <- sum(sort(out_preds2)[39:75]*
                  analysisdat2$cc.wts[order(out_preds2)[39:75]])
pi2_sum3 <- sum(sort(out_preds2)[76:112]*
                  analysisdat2$cc.wts[order(out_preds2)[76:112]])
pi2_sum4 <- sum(sort(out_preds2)[113:150]*
                  analysisdat2$cc.wts[order(out_preds2)[113:150]])

y2_sum1 <- sum(analysisdat2$Y[order(out_preds2)[1:38]]*
                 analysisdat2$cc.wts[order(out_preds2)[1:38]])
y2_sum2 <- sum(analysisdat2$Y[order(out_preds2)[39:75]]*
                 analysisdat2$cc.wts[order(out_preds2)[39:75]])
y2_sum3 <- sum(analysisdat2$Y[order(out_preds2)[76:112]]*
                 analysisdat2$cc.wts[order(out_preds2)[76:112]])
y2_sum4 <- sum(analysisdat2$Y[order(out_preds2)[113:148]]*
                 analysisdat2$cc.wts[order(out_preds2)[113:148]])

ts2 <- (y2_sum1 - pi2_sum1)^2 / (pi2_sum1*(1 - pi2_sum1/38)) +
  (y2_sum2 - pi2_sum2)^2 / (pi2_sum2*(1 - pi2_sum2/37)) +
  (y2_sum3 - pi2_sum3)^2 / (pi2_sum3*(1 - pi2_sum3/37)) +
  (y2_sum4 - pi2_sum4)^2 / (pi2_sum4*(1 - pi2_sum4/38))

pchisq(ts2, 2, lower.tail = FALSE)


#########################################################################
# Sensitivity analysis for RII removing two outliers based on diagnostics
#########################################################################
analysisdat2_sens <- analysisdat2 %>%
  filter(ID %in% analysisdat2$ID[prop_mod2$residuals > -3])

# Initialize DR estimates and confidence intervals
dr_ests2 <- array(NA, dim = c(4, length(m2vals)))
dr_lower2 <- array(NA, dim = c(4, length(m2vals)))
dr_upper2 <- array(NA, dim = c(4, length(m2vals)))

# Fit DR model with no measurement error
dr.res2.0 <- fit.dr(data = analysisdat2_sens, args = dr.args2, a = m2vals)

# Grab results corresponding to counterfactual curve
EYa2.ind <- grepl("EYa", names(dr.res2.0$est))
EYa2 <- dr.res2.0$est[EYa2.ind]
cov.EYa2 <- dr.res2.0$bc.var[EYa2.ind, EYa2.ind]

# Delta method to log
grad.log2 <- diag(1 / EYa2)
se.log.EYa2 <- sqrt(diag(grad.log2 %*% cov.EYa2 %*% t(grad.log2)))

# Add results with no measurement error
dr_ests2[1, ] <- EYa2
dr_lower2[1, ] <- exp(log(EYa2) - qnorm(0.975) * se.log.EYa2)
dr_upper2[1, ] <- exp(log(EYa2) + qnorm(0.975) * se.log.EYa2)

# Results for increasing measurement error
for (me in 2:4) {

  sigma_me2 <- var(analysisdat2_sens$A)*me_list2[me]

  # DR estimator
  dr.res2 <- fit.dr.mccs(data = analysisdat2_sens, a = m2vals,
                         cov.e = sigma_me2, B = 100, mc.seed = 12345,
                         return.var = TRUE, args = dr.args2,
                         start = dr_ests2[me - 1, 1:8])

  # Grab results corresponding to counterfactual curve
  EYa2.ind <- grepl("EYa", names(dr.res2$est))
  EYa2 <- dr.res2$est[EYa2.ind]
  cov.EYa2 <- dr.res2$bc.var[EYa2.ind, EYa2.ind]

  # Delta method to log
  grad.log2 <- diag(1 / EYa2)
  se.log.EYa2 <- sqrt(diag(grad.log2 %*% cov.EYa2 %*% t(grad.log2)))
  dr_ests2[me, ] <- EYa2
  dr_lower2[me, ] <- exp(log(EYa2) - qnorm(0.975) * se.log.EYa2)
  dr_upper2[me, ] <- exp(log(EYa2) + qnorm(0.975) * se.log.EYa2)

}

# Plot sensitivity analysis results
latdat_sens <- data.frame(vals = c(rep(m1vals, 4),
                              rep(m2vals, 4)),
                     Risk = c(dr_ests1[1, ], dr_ests1[2, ],
                              dr_ests1[3, ], dr_ests1[4, ],
                              dr_ests2[1, ], dr_ests2[2, ],
                              dr_ests2[3, ], dr_ests2[4, ]),
                     Risk_low = c(dr_lower1[1, ], dr_lower1[2, ],
                                  dr_lower1[3, ], dr_lower1[4, ],
                                  dr_lower2[1, ], dr_lower2[2, ],
                                  dr_lower2[3, ], dr_lower2[4, ]),
                     Risk_upp = c(dr_upper1[1, ], dr_upper1[2, ],
                                  dr_upper1[3, ], dr_upper1[4, ],
                                  dr_upper2[1, ], dr_upper2[2, ],
                                  dr_upper2[3, ], dr_upper2[4, ]),
                     ME = c(rep(round(me_list, 3), each = length(m1vals)),
                            rep(round(me_list2, 3), each = length(m2vals))),
                     Exposure = c(rep("ADCP", 4*length(m1vals)),
                                  rep("RII", 4*length(m2vals))))
latdat_sens$Risk_upp[latdat_sens$Risk_upp > 1] <- 1

# New facet label names for Exposure variable
ggplot(latdat_sens, aes(x = vals, y = Risk, ymin = Risk_low, ymax = Risk_upp)) +
  geom_line() +
  facet_grid(ME ~ Exposure, scales = "free",
             labeller = label_bquote(sigma[me]^2 == .(ME) * sigma^2)) +
  #labeller = labeller(ME = me.labs, Exposure = exp.labs)) +
  geom_ribbon(alpha = 0.3) +
  xlab("Exposure values") + ylab("HIV risk") +
  ylim(c(0, 1)) +
  theme_bw()

# Updated diagnostics for sensitivity analysis
prop_mod2_sens <- lm(A ~ L1 + L2 + L3 + L4 + L5 + L6, data = analysisdat2_sens)
out_mod2_sens <- fit.glm(analysisdat2_sens, dr.args2)
out_preds2_sens <- exp(out_mod2_sens$est %*% t(cbind(1, analysisdat2_sens[, 2:8])))

# Propensity model diagnostics
par(mfrow = c(2, 2))
plot(prop_mod2_sens)

# Outcome model diagnostics (manual Chi-square GoF)
# RII outcome model
pi2_sum1 <- sum(sort(out_preds2_sens)[1:37]*
                  analysisdat2_sens$cc.wts[order(out_preds2_sens)[1:37]])
pi2_sum2 <- sum(sort(out_preds2_sens)[38:74]*
                  analysisdat2_sens$cc.wts[order(out_preds2_sens)[38:74]])
pi2_sum3 <- sum(sort(out_preds2_sens)[75:111]*
                  analysisdat2_sens$cc.wts[order(out_preds2_sens)[75:111]])
pi2_sum4 <- sum(sort(out_preds2_sens)[112:148]*
                  analysisdat2_sens$cc.wts[order(out_preds2_sens)[112:148]])

y2_sum1 <- sum(analysisdat2_sens$Y[order(out_preds2_sens)[1:37]]*
                 analysisdat2_sens$cc.wts[order(out_preds2_sens)[1:37]])
y2_sum2 <- sum(analysisdat2_sens$Y[order(out_preds2_sens)[38:74]]*
                 analysisdat2_sens$cc.wts[order(out_preds2_sens)[38:74]])
y2_sum3 <- sum(analysisdat2_sens$Y[order(out_preds2_sens)[75:111]]*
                 analysisdat2_sens$cc.wts[order(out_preds2_sens)[75:111]])
y2_sum4 <- sum(analysisdat2_sens$Y[order(out_preds2_sens)[112:148]]*
                 analysisdat2_sens$cc.wts[order(out_preds2_sens)[112:148]])

ts2 <- (y2_sum1 - pi2_sum1)^2 / (pi2_sum1*(1 - pi2_sum1/37)) +
  (y2_sum2 - pi2_sum2)^2 / (pi2_sum2*(1 - pi2_sum2/37)) +
  (y2_sum3 - pi2_sum3)^2 / (pi2_sum3*(1 - pi2_sum3/37)) +
  (y2_sum4 - pi2_sum4)^2 / (pi2_sum4*(1 - pi2_sum4/37))

pchisq(ts2, 2, lower.tail = FALSE)
