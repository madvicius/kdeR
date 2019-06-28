

.shuffle_rows <- function(mat) {
    mat[sample(nrow(mat)), ]
}

#' Amostra uma mistura de normais univariada
#' @export rnormmix
#' @param n número de observacões
#' @param mu vetor das médias
#' @param sigma vetor da variâncias
#' @param prop vetor das probabilidades
#' @examples
#'  x <- rnormmix(100, c(-5, 5))
#' xy <- kde(x, h='STE')
#' plot(xy)
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
#' @examples
#' sim <- rnormmix(105, mu = c(-1.75, 0.6, 3.0),
#' sigma = c(0.75, 0.5, 0.75), prob = c(3/6, 1/6, 2/6))
#'
#' dsim <- normmix(sort(sim1), mu = c(-1.75, 0.6, 3.0),
#' sigma = c(0.75, 0.5, 0.75), prob = c(3/6, 1/6, 2/6))
#'
#' plot(sort(sim),dsim, type = 'l')
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
#' @examples
#' xy <- rbnormmix(100, list(c(-5, -5), c(5, 5)))
#'
#' bmix_bkde <- bkde(xy, h='STE')
#'
#' persp(bmix_bkde,theta=-120,phi=30)
#' plot(bmix_bkde,pch=16)
#' contour(bmix_bkde,add=TRUE)
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
