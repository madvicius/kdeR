

#' Estima a distribuição univariada de uma amostra por KDE
kde <- function(
    sp, h = c('NS', 'STE'),
    kernel = dnorm, ngrid = 10000
) {
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

#' Plota a densidade estimada por KDE
plot.kde <- function(kde, xlab='x', ylab='density',...) {
    plot(kde$x, kde$y, type='l', xlab=xlab, ylab=ylab,... )
}

#' Estima a distribuição bivariada de uma amostra por KDE
bkde <- function(
    sp, h = c('NS', 'STE'),
    kernel = dmvnorm, ngrid = 100
) {
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

    
    
    