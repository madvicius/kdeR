

.shuffle_rows <- function(mat) {
    mat[sample(nrow(mat)), ]
}

#' Amostra uma mistura de normais univariada
#' @export rnormmix
#' @param n número de observacões
#' @param mu vetor das médias
#' @param sigma vetor da variâncias
#' @param prop vetor das probabilidades
rnormmix <- function(
    n,
    mu,
    sigma = rep(1, length(mu)),
    prob = rep(1/length(mu), length(mu))
) {

  require(purrr)
  require(dplyr)
    ni <- as.vector(rmultinom(1, n, prob))

    list(ni, mu, sigma) %>%
        pmap(~ rnorm(..1, ..2, ..3)) %>%
        simplify() %>%
        sample()
}

#' Densidade da mistura de normais univariada
#' @export dnormmix
#' @param x vetor de quantis
#' @param mu vetor das médias
#' @param sigma vetor da variâncias
#' @param prop vetor das probabilidades
dnormmix <- function(
    x,
    mu,
    sigma = rep(1, length(mu)),
    prob = rep(1/length(mu), length(mu))
) {

  require(purrr)
  require(dplyr)
    x %>% map_dbl(~ sum(prob * dnorm(.x, mu, sigma)))
}

#' Amostra uma mistura de normais bivariadas
#' @export rbnormmix
#' @param n número de observacões
#' @param mu vetor das médias
#' @param sigma vetor da variâncias
#' @param prop vetor das probabilidades
rbnormmix <- function(
    n,
    mu,
    sigma = rep(list(diag(length(mu))), length(mu)),
    prob = rep(1/length(mu), length(mu))
) {
  require(mvtnorm)
  require(purrr)
  require(dplyr)
    ni <- as.vector(rmultinom(1, n, prob/sum(prob)))

    list(ni, mu, sigma) %>%
        pmap(~ rmvnorm(..1, ..2, ..3)) %>%
        do.call(rbind, .) %>%
        .shuffle_rows()
}
