


#' Estima a distribuicão univariada de uma amostra por KDE
#' @export kde
#' @param sp vetor numérico de qual o kde irá estimar a densidade
#' @param h método utilizado para o cálculo de h. A funcao possui duas opcões
#' , "STE" (solve the equation) e "NS" (normal scale). default = STE
#' @param kernel kernel utilizado na estimacão. default = dnorm
#' @param ngrid numero de pontos na malha. default = 10000
#' @examples
#' x <- rnormmix(100, c(-5, 5))
#' xy <- kde(x, h='STE')
#'
#' plot(xy)
kde <- function(
    sp, h = c('NS', 'STE'),
    kernel = dnorm, ngrid = 10000
) {
  require(purrr)
  require(dplyr)
    n <- length(sp)
    if (is.character(h)) { h <- switch( h[1],
        `NS` = bw.nrd(sp),
        `STE` = bw.SJ(sp)
    )}

    fhat <- function(x) {
        x %>%
            map(~ (.x - sp) / h) %>%
            map(kernel) %>%
            map(sum) %>%
            map(~ .x / (n*h)) %>%
            simplify()
    }
    spy <- fhat(sp)

    x <- seq(
        min(sp) - 5*h,
        max(sp) + 5*h,
        length.out = ngrid
    )
    y <- fhat(x)

    structure(list(
        x = x, y = y, fhat = fhat,
        sp = sp, spy = spy, h = h,
        kernel = kernel, n = n, ngrid = ngrid),
        class = 'kde'
    )
}

#' Plota a densidade estimada pelo KDE univariado
#' @export plot.kde
#' @param kde objeto da classe kde
#' @param ... você pode utilizar outros argumentos da funcão plot.default
plot.kde <- function(kde, xlab='x', ylab='density',...) {
    plot(kde$x, kde$y, type='l', xlab=xlab, ylab=ylab,... )
}

#' Estima a distribuicão bivariada de uma amostra por KDE
#' @export bkde
#' @param sp vetor numérico de qual o bkde irá estimar a densidade
#' @param h método utilizado para o cálculo de h. A funcao possui duas opcões
#' , "STE" (solve the equation) e "NS" (normal scale). default = STE
#' @param kernel kernel utilizado na estimacão. default = dmvnorm
#' @param ngrid numero de pontos na malha. default = 100
#' @examples
#' xy <- rbnormmix(100, list(c(-5, -5), c(5, 5)))
#'
#' bmix_bkde <- bkde(xy, h='STE')
#'
#' plot(bmix_bkde,pch=16)
#' contour(bmix_bkde,add=TRUE)
bkde <- function(
    sp, h = c('NS', 'STE'),
    kernel = dmvnorm, ngrid = 100
) {
  require(purrr)
  require(dplyr)
    n <- nrow(sp)
    if (is.character(h)) {
        h1 <- switch( h[1],
            `NS` = bw.nrd(sp[,1]),
            `STE` = bw.SJ(sp[,1])
        )
        h2 <- switch( h[1],
            `NS` = bw.nrd(sp[,2]),
            `STE` = bw.SJ(sp[,2])
        )
    } else {
        h1 <- h[1]
        h2 <- h[2]
    }

    fhat <- function(x, y) {
        list(x, y) %>%
            pmap(~ cbind(
                (.x - sp[,1]) / h1,
                (.y - sp[,2]) / h2
            )) %>%
            map(kernel) %>%
            map(sum) %>%
            map(~ .x / (n*h1*h2)) %>%
            simplify()
    }
    spz <- fhat(sp[,1], sp[,2])

    x <- seq(
        min(sp[,1]) - 5*h1,
        max(sp[,1]) + 5*h1,
        length.out = ngrid
    )
    y <- seq(
        min(sp[,2]) - 5*h2,
        max(sp[,2]) + 5*h2,
        length.out = ngrid
    )
    xy <- expand.grid(x, y)
    z <- matrix(fhat(xy[,1], xy[,2]), nrow = ngrid)

    structure(list(
        x = x, y = y, z = z, fhat = fhat,
        sp = sp, spz = spz, h1 = h1, h2 = h2,
        kernel = kernel, n = n, ngrid = ngrid),
        class = 'bkde'
    )
}


#' Plota os dados de BKDE
#' @export plot.bkde
#' @param bkde objeto da classe bkde
#' @param ... você pode utilizar outros argumentos da funcão plot.default
#' @examples
#'  xy <- rbnormmix(100, list(c(-5, -5), c(5, 5)))
#'
#' bmix_bkde <- bkde(xy, h='STE')
#'
#' plot(bmix_bkde,pch=16)
#' contour(bmix_bkde,add=TRUE)
plot.bkde <- function(bkde,
                      xlab= 'x',
                      ylab= 'y',
                      ...){
    plot(x=bkde$sp[,1],
         y=bkde$sp[,2],
         xlab = xlab,
         ylab= ylab,
         ...)
}
#' Plota um gráfico de contorno de BKDE
#' @export contour.bkde
#' @param bkde objeto da classe bkde
#' @param quantiles defini os niveis do gráfico de contorno
#' @parm ... você pode utilizar os outros parametros da função contour.default
#' @examples
#' xy <- rbnormmix(100, list(c(-5, -5), c(5, 5)))
#'
#' bmix_bkde <- bkde(xy, h='STE')
#'
#' plot(bmix_bkde,pch=16)
#' contour(bmix_bkde,add=TRUE)

contour.bkde <- function(bkde,
                         xlab= 'x',
                         ylab= 'y',
                         quantiles=c(0,.25,.5),
                         drawlabels=FALSE,
                         ...){

    quantiles <- quantile(bkde$spz, quantiles)
    quantiles[1] <- quantiles[1]/2
        contour(bkde$z,
            x= bkde$x,
            y= bkde$y,
            levels = quantiles,
            drawlabels = drawlabels,
            ...)
}

#' Plota a densidade estimada por BKDE
#' @export persp.bkde
#' @param bkde objeto da classe bkde
#' @parm ... você pode utilizar os outros parametros da função persp.default
#' @examples
#' xy <- rbnormmix(100, list(c(-5, -5), c(5, 5)))
#'
#' bmix_bkde <- bkde(xy, h='STE')
#'
#' persp(bmix_bkde,theta=-120,phi=30)

persp.bkde <- function(bkde,
                      xlab='x',
                      ylab='y',
                      zlab = 'density',
                      ticktype = 'detailed',
                      ...){
    persp(z= bkde$z,
          x= bkde$x,
          y= bkde$y,
          ticktype = ticktype,
          zlab = zlab,
          xlab= xlab,
          ylab= ylab,
          ...
    )
}



