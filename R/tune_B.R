#' Tune MCCS estimating equation to select B value
#'
#'
#' @export
tune.B <- function(get.psi, data, cov.e, BB,
                   args = list(), mc.seed = 123) {

  ## for troubleshooting
  #get.psi <- get.psi.glm; data <- datstar; BB <- c(5, 10, 20)

  ## unpack arguments
  list2env(args, envir = environment())

  ## loop through B values and evaluate estimating equation
  tab <- t(vapply(
    X = BB,
    FUN.VALUE = c(0, g),
    FUN = function(B) {

      ## create MCCS GLM estimating function
      get.psi.mccs <- make.mccs(
        get.psi = get.psi, data = data, args = args,
        cov.e = cov.e, B = B, mc.seed = mc.seed)

      c(B = B,
        psi = get.psi.mccs(x = g))
    }
  )) %>%
    as.data.frame()

  ## plot results
  plt <- tab %>%
    pivot_longer(cols = !B,
                 names_to = "Component",
                 values_to = "Psi") %>%
    ggplot(aes(x = B,
               y = Psi,
               color = Component)) +
    geom_line()

  return(list(data = tab,
              plot = plt))

}
