#' Locally change-point detection using the Mann-Whitney test
#'
#' Locally change-point detection after accommodating the Mann-Whitney test in the LACPD procedure
#'
#'
#'
#' @param x a numeric vector
#' @param m number of times to sub-sample
#' @param k single number or numeric vector proportional to the number of points on each side of the target point. See details
#' @param blow fraction of observations (0-1) at the beginning of the time-series not considered for change detection
#' @param bup similar to \code{blow}, but for the end of the time series. Default is 1-\code{blow}
#' @param adjust if \code{TRUE}, p-value will be adjusted by methods in \link[stats]{p.adjust}
#' @param history if \code{TRUE}, it maintains the step-wise results when \code{k} is a vector
#' @param eps an integer that will be used as a stopping rule. Once the difference between three different sliding windows in a row is less than eps, the algorithm stops.
#' @param alpha significance level
#' @param double logical. if TRUE it adjust p-values two times.
#' @param ... arguments passed to \link[stats]{p.adjust}
#'
#'
#' @details
#' This technique accommodates the \link[stats]{wilcox.test} in the LACPD procedure of Moradi et al. (2020) to look for potential change-points in the numerical vector \code{x}.
#'
#' In the tails of \code{x}, since there are not enough data points before/after the target points, the method uses sub-sampling, wherein it moderates the effect of sub-sampling by iteration.
#'
#' Assume the length of \code{x} is n. The argument \code{k} is used to set the number of data in the sides of the target point when looking for change-points. n/\code{k}  is the number of points on each side we consider. For instance, if n=300 and k=10, this means we consider 30 observations before and 30 after when locally detecting changes.
#'
#' If \code{adjust=TRUE}, methods susch as "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" can be passed through, e.g. method="BY".
#'
#' If \code{k} is a vector of numbers, then the function returns a result based on adaptive sliding windows which is an average result obtained from different windows.
#'
#'
#'
#' @return
#' adaptive window: the automatic selected adaptive window.
#'
#' change-point: the index of the most probable change point in the time series
#'
#' statistics: the z statistics of LACPD
#'
#' magnitude: the magnitude of change
#'
#' p.value: the corresponding p.value
#'
#'
#' Attributes:
#'
#' attr(,"zs"): retrieves the obtained z statistics at the time-periods the function has looked for potential changes
#'
#' attr(,"ps"): retrieves the obtained p-values at the time-periods the function has looked for potential changes
#'
#' attr(,"mags"): retrieves the magnitude of change at the time-periods the function has looked for potential changes
#'
#'
#' if \code{history=TRUE} and \code{k} is a numeric vector, then the following attributes can also be retrieved.
#'
#' attr(,"history"): a dataframe containing the results (cp, magnitude, Z, and p.value) based on different adaptive sliding windows which are permutations of the given \code{k}
#'
#' attr(,"allzs"): a list which retrieves the obtained z statistics at the time-periods the function has looked for potential changes based on different adaptive sliding windows which are permutations of the given \code{k}
#'
#' attr(,"allmags"): a list which retrieves the magnitude of change at the time-periods the function has looked for potential changes based on different adaptive sliding windows which are permutations of the given \code{k}
#'
#' attr(,"allps"): a list which retrieves the obtained p-values at the time-periods the function has looked for potential changes based on different adaptive sliding windows which are permutations of the given \code{k}
#'
#' attr(,"zmat"): a matrix of z statistics for each individual sliding window
#'
#' attr(,"pmat"): a matrix of p-values for each individual sliding window
#'
#' attr(,"mmat"): a matrix of magnitudes for each individual sliding window
#'
#'
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#'
#' @references
#' Moradi, M., Montesino-SanMartin, M., Ugarte, M. D., and Militino, A. F. (2020). Locally adaptive change-point detection with applications to remote sensing and land use changes.
#'
#' @seealso
#' \link[stats]{wilcox.test}, \link[stats]{p.adjust}
#'
#' @examples
#' x <- rnorm(50)
#' Z <- lacpd_mw(x,m=10,k=3)
#' plot(Z$s,attr(Z,"zs"),type="l",ylab = "Z",xlab = "time")
#' plot(Z$s,attr(Z,"mags"),type="l",ylab = "Z",xlab = "time")


#' @import stats
#' @export
lacpd_mw <- function(x,m=1,k=c(2:5),blow=0.1,bup=(1-blow),
                     adjust=FALSE,history=FALSE,eps=1,alpha=0.05,double=FALSE,...){
  n <- length(x)

  a <- round(blow*n)
  b <- round(bup*n)
  s <- seq(a,b,by=1)


  if(anyNA(x)){

    p.val <- NA
    cp <- NA
    mag <- NA
    z <- NA


    out.final <- list(s=s,
                      cp=cp,
                      z=z,
                      mag=mag,
                      p=p.val)
    class(out.final) <- c("list","lacpd.mw")
    return(out.final)
  }


  ############### if k contains several values

  if(length(k)>1){

    out.list.z <- list()
    out.list.p <- list()
    out.list.mag <- list()

    if(history){ hist1 <- list() }

    for (i in 1:length(k)){

      out.pre <- lacpd_mw(x=x,m=m,k=k[i],blow=blow,bup=(1-blow),adjust=adjust,...)


      out.list.z[[i]] <- attr(out.pre,"zs")
      out.list.p[[i]] <- attr(out.pre,"ps")
      out.list.mag[[i]] <- attr(out.pre,"mags")

      if(history){hist1[[i]] <- c(out.pre$cp,out.pre$mag,out.pre$z,out.pre$p)}

      if(i>3){

        if(double){
          pcheck <- p.adjust(colMeans(do.call(rbind,out.list.p[1:i])),...)
        }else{
          pcheck <- colMeans(do.call(rbind,out.list.p[1:i]))
        }

        if(min(pcheck)>alpha){

          win <- paste(k[1],"-",k[i-1])

          if(history){
            hist1 <- do.call(rbind,hist1)
            hist1 <- hist1[c(1:(i-1)),]
          }

          out.list.z <- do.call(rbind,out.list.z)
          out.list.z <- out.list.z[c(1:(i-1)),]
          zs <- colMeans(out.list.z)

          out.list.p <- do.call(rbind,out.list.p)
          out.list.p <- out.list.p[c(1:(i-1)),]
          ps <- colMeans(out.list.p)

          out.list.mag <- do.call(rbind,out.list.mag)
          out.list.mag <- out.list.mag[c(1:(i-1)),]
          mags <- colMeans(out.list.mag)

          break()

        }else if(
          (max(abs(s[which.min(colMeans(do.call(rbind,out.list.p[1:i])))]
                   -
                   s[which.min(colMeans(do.call(rbind,out.list.p[1:(i-1)])))])
               ,
               abs(s[which.min(colMeans(do.call(rbind,out.list.p[1:i])))]
                   -
                   s[which.min(colMeans(do.call(rbind,out.list.p[1:(i-2)])))])) < eps)
        ){

          win <- paste(k[1],"-",k[(i-1)])

          if(history){
            hist1 <- do.call(rbind,hist1)
            hist1 <- hist1[c(1:(i-1)),]
          }


          out.list.z <- do.call(rbind,out.list.z)
          out.list.z <- out.list.z[c(1:(i-1)),]
          zs <- colMeans(out.list.z)

          out.list.p <- do.call(rbind,out.list.p)
          out.list.p <- out.list.p[c(1:(i-1)),]
          ps <- colMeans(out.list.p)

          out.list.mag <- do.call(rbind,out.list.mag)
          out.list.mag <- out.list.mag[c(1:(i-1)),]
          mags <- colMeans(out.list.mag)

          break()

        }
      }

      if(i==length(k)){

        win <- paste(k[1],"-",k[(i)])

        if(history){hist1 <- do.call(rbind,hist1)}
        # hist1 <- hist1[c(1:(i-2)),]

        out.list.z <- do.call(rbind,out.list.z)
        # out.list.z <- out.list.z[c(1:(i-2)),]
        zs <- colMeans(out.list.z)

        out.list.p <- do.call(rbind,out.list.p)
        # out.list.p <- out.list.p[c(1:(i-2)),]
        ps <- colMeans(out.list.p)

        out.list.mag <- do.call(rbind,out.list.mag)
        # out.list.mag <- out.list.mag[c(1:(i-2)),]
        mags <- colMeans(out.list.mag)
      }
    }
    # win <- paste(k[1],"-",k[(i)])
    # hist1 <- do.call(rbind,hist1)
    # hist1 <- hist1[-nrow(hist1),]
    #
    # out.list.z <- do.call(rbind,out.list.z)
    # # out.list.z <- out.list.z[-nrow(out.list.z),]
    # zs <- colMeans(out.list.z)
    # #
    # out.list.p <- do.call(rbind,out.list.p)
    # # out.list.p <- out.list.p[-nrow(out.list.p),]
    # ps <- colMeans(out.list.p)
    # #
    # out.list.mag <- do.call(rbind,out.list.mag)
    # # out.list.mag <- out.list.mag[-nrow(out.list.mag),]
    # mags <- colMeans(out.list.mag)

    if(history){

      out.hist.zs <- list()
      out.hist.ps <- list()
      out.hist.mags <- list()
      out.hist <- list()

      rn1 <- c()
      for (i in 1:nrow(hist1)) {
        rn1[i] <- paste(k[i])
      }
      colnames(hist1) <- c("change-point","magnitude","Z","p.value")
      rownames(hist1) <- rn1

      for (h in 1:(nrow(out.list.z)-1)) {

        allzs <- list()
        allmags <- list()
        allps <- list()

        cphist <- c()
        maghist <- c()
        zhist <- c()
        phist <- c()

        for (i in (h+1):nrow(out.list.z)){

          allzs[[i-h]] <- colMeans(out.list.z[h:i,])
          allps[[i-h]] <- colMeans(out.list.p[h:i,])
          allmags[[i-h]] <- colMeans(out.list.mag[h:i,])

          # if(adjust){
          #   allps[[i-h]] <- p.adjust(allps[[i-h]],...)
          # }else{
          #   allps[[i-h]] <- allps[[i-h]]
          # }

          pcod <- which.min(allps[[i-h]])
          cphist[i-h] <- s[pcod]
          maghist[i-h] <- allmags[[i-h]][pcod]
          zhist[i-h] <- allzs[[i-h]][pcod]

          phist[i-h] <- min(allps[[i-h]])

        }

        hist <- cbind(cphist,maghist,zhist,phist)

        rn <- c()
        for (i in 1:nrow(hist)) {
          rn[i] <- paste(k[h],"-",k[i+h])
        }
        colnames(hist) <- c("change-point","magnitude","Z","p.value")
        rownames(hist) <- rn

        names(allps) <- names(allmags) <- names(allzs) <- rn


        out.hist[[h]] <- hist
        out.hist.zs[[h]] <- allzs
        out.hist.mags[[h]] <- allmags
        out.hist.ps[[h]] <- allps

      }



      out.hist.zs <- unlist(out.hist.zs, recursive=FALSE)
      out.hist.ps <- unlist(out.hist.ps, recursive=FALSE)
      out.hist.mags <- unlist(out.hist.mags, recursive=FALSE)
      out.hist <- do.call(rbind,out.hist)

    }

    if(history){out.hist <- rbind(hist1,out.hist)}

    if(double){
      ps <- p.adjust(ps,...)
    }

    cp <- s[which.min(ps)]
    mag <- mags[which.min(ps)]
    z <- zs[which.min(ps)]

    p <- min(ps)


    out.final <- list(window=win,
                      s=s,
                      cp=cp,
                      z=z,
                      mag=mag,
                      p=p)

    class(out.final) <- c("list","lacpd.mw")
    attr(out.final,"zs") <- zs
    attr(out.final,"mags") <- mags
    attr(out.final,"ps") <- ps

    attr(out.final,"zmat") <- out.list.z
    attr(out.final,"pmat") <- out.list.p
    attr(out.final,"mmat") <- out.list.mag


    if(history){
      attr(out.final,"history") <- out.hist
      attr(out.final,"allzs") <- out.hist.zs
      attr(out.final,"allmags") <- out.hist.mags
      attr(out.final,"allps") <- out.hist.ps
    }

    return(out.final)
  }



  low <- ceiling(n/k)
  up <- floor(n-(n/k))
  if(low < a) low <- a #stop(" n/k should be larger than blow*n where n is the length of x ")
  if(up > b)  up <- b

  low1 <- which(s==low)
  up1 <- which(s==up)



  ############### if k>2, then subsampling is done on the first and last part of x

  if(k==2){

    ################################################ lower part

    out.low <- lapply(X=1:m, function(i){

      zs <- numeric()
      mag <- numeric()
      ps <- numeric()

      for (j in 1:(low1)){

        z <- c(sample(x[1:(s[j]-1)],size = (n-s[j]+1),replace = T),
               x,
               sample(x[(s[j]+1):n],size =(s[j]),replace = T ))

        znew <- z[ceiling(((n*(k-1))/k)+1):(n+(n/k)+1)]


        znew <- znew[-(floor(length(znew)/2)+1)]
        zlen <- length(znew)
        mag[j] <- abs(mean(znew[((low+1)+1):zlen]) - mean(znew[1:(low)]))

        mk <- wilcox.test(znew[1:low],znew[(low+1):zlen])
        zs[j] <- as.numeric(mk$statistic)
        ps[j] <- as.numeric(mk$p.value)

      }

      return(list(zs=zs,mag=mag,ps=ps))

    })

    out.z.low <- lapply(X=1:m, function(i){
      out.low[[i]]$zs
    })
    zs.low <- rowMeans(do.call(cbind,out.z.low))

    out.p.low <- lapply(X=1:m, function(i){
      out.low[[i]]$ps
    })
    ps.low <- rowMeans(do.call(cbind,out.p.low))

    out.mag.low <- lapply(X=1:m, function(i){
      out.low[[i]]$mag
    })
    mag.low <- rowMeans(do.call(cbind,out.mag.low))

    ################################################## upper part
    out.up <- lapply(X=1:m, function(i){
      zs <- numeric()
      mag <- numeric()
      ps <- numeric()

      for (j in (up1+1):length(s)){

        z <- c(sample(x[1:(s[j]-1)],size = (n-s[j]+1),replace = T),
               x,
               sample(x[(s[j]+1):n],size =(s[j]),replace = T ))

        znew <- z[ceiling(((n*(k-1))/k)+1):(n+(n/k)+1)]


        znew <- znew[-(floor(length(znew)/2)+1)]
        zlen <- length(znew)
        mag[j] <- abs(mean(znew[(up+1):zlen]) - mean(znew[1:up]))


        mk <- wilcox.test(znew[1:up],znew[(up+1):zlen])
        zs[j] <- as.numeric(mk$statistic)
        ps[j] <- as.numeric(mk$p.value)

      }

      return(list(zs=zs[(up1+1):length(s)],mag=mag[(up1+1):length(s)],
                  ps=ps[(up1+1):length(s)]))
    })

    out.z.up <- lapply(X=1:m, function(i){
      out.up[[i]]$zs
    })
    zs.up <- rowMeans(do.call(cbind,out.z.up))

    out.p.up <- lapply(X=1:m, function(i){
      out.up[[i]]$ps
    })
    ps.up <- rowMeans(do.call(cbind,out.p.up))

    out.mag.up <- lapply(X=1:m, function(i){
      out.up[[i]]$mag
    })
    mag.up <- rowMeans(do.call(cbind,out.mag.up))


    zs <- c(zs.low,zs.up)
    ps <- c(ps.low,ps.up)
    mags <- c(mag.low,mag.up)

  }
  else{ ########## k>2
    ################################################ lower part
    out.low <- lapply(X=1:m, function(i){
      zs <- numeric()
      mag <- numeric()
      ps <- numeric()

      for (j in 1:(low1)){

        z <- c(sample(x[1:(s[j]-1)],size = (n-s[j]+1),replace = T),
               x,
               sample(x[(s[j]+1):n],size =(s[j]),replace = T ))
        znew <- z[ceiling(((n*(k-1))/k)+1):(n+(n/k)+1)]


        znew <- znew[-(floor(length(znew)/2)+1)]
        zlen <- length(znew)
        mag[j] <- abs(mean(znew[((zlen/2)+1):zlen]) - mean(znew[1:(zlen/2)]))

        midz <- zlen/2
        mk <- wilcox.test(znew[1:midz],znew[(midz+1):zlen])
        zs[j] <- as.numeric(mk$statistic)
        ps[j] <- as.numeric(mk$p.value)

      }

      return(list(zs=zs,mag=mag,ps=ps))

    })

    out.z.low <- lapply(X=1:m, function(i){
      out.low[[i]]$zs
    })
    zs.low <- rowMeans(do.call(cbind,out.z.low))

    out.p.low <- lapply(X=1:m, function(i){
      out.low[[i]]$ps
    })
    ps.low <- rowMeans(do.call(cbind,out.p.low))

    out.mag.low <- lapply(X=1:m, function(i){
      out.low[[i]]$mag
    })
    mag.low <- rowMeans(do.call(cbind,out.mag.low))

    ################################################ middle part
    zs.middle <-  numeric()
    ps.middle <- numeric()
    mag_m <- numeric()
    for (j in (low1+1):(up1-1)){

      z <- c(sample(x[1:(s[j]-1)],size = (n-s[j]+1),replace = T),
             x,
             sample(x[(s[j]+1):n],size =(s[j]),replace = T ))

      znew <- z[ceiling(((n*(k-1))/k)+1):(n+(n/k)+1)]


      znew <- znew[-(floor(length(znew)/2)+1)]
      zlen <- length(znew)
      mag_m[j] <- abs(mean(znew[((zlen/2)+1):zlen]) - mean(znew[1:(zlen/2)]))

      midz <- zlen/2
      mk <-  wilcox.test(znew[1:midz],znew[(midz+1):zlen])
      zs.middle[j] <- as.numeric(mk$statistic)
      ps.middle[j] <- as.numeric(mk$p.value)

    }
    zs.middle <- zs.middle[(low1+1):(up1-1)]
    ps.middle <- ps.middle[(low1+1):(up1-1)]
    mag.middle <- mag_m[(low1+1):(up1-1)]
    ################################################ upper part
    out.up <- lapply(X=1:m, function(i){

      zs <- numeric()
      ps <- numeric()
      mag <- numeric()

      for (j in up1:length(s)){
        z <- c(sample(x[1:(s[j]-1)],size = (n-s[j]+1),replace = T),
               x,sample(x[(s[j]+1):n],size =(s[j]),replace = T ))

        znew <- z[ceiling(((n*(k-1))/k)+1):(n+(n/k)+1)]


        znew <- znew[-(floor(length(znew)/2)+1)]
        zlen <- length(znew)
        mag[j] <- abs(mean(znew[((zlen/2)+1):zlen]) - mean(znew[1:(zlen/2)]))

        midz <- zlen/2
        mk <-  wilcox.test(znew[1:midz],znew[(midz+1):zlen])
        zs[j] <- as.numeric(mk$statistic)
        ps[j] <- as.numeric(mk$p.value)
      }
      return(list(zs=zs[up1:length(s)],mag=mag[up1:length(s)],
                  ps=ps[up1:length(s)]))
    })

    out.z.up <- lapply(X=1:m, function(i){
      out.up[[i]]$zs
    })
    zs.up <- rowMeans(do.call(cbind,out.z.up))

    out.p.up <- lapply(X=1:m, function(i){
      out.up[[i]]$ps
    })
    ps.up <- rowMeans(do.call(cbind,out.p.up))

    out.mag.up <- lapply(X=1:m, function(i){
      out.up[[i]]$mag
    })
    mag.up <- rowMeans(do.call(cbind,out.mag.up))

    zs <- c(zs.low,zs.middle,zs.up)
    ps <- c(ps.low,ps.middle,ps.up)
    mags <- c(mag.low,mag.middle,mag.up)

  }

  if(adjust){
    ps <- p.adjust(ps,...)
  }else{
    ps <- ps
  }

  p.val <- min(ps)
  cp <- s[which.min(ps)]
  mag <- mags[which.min(ps)]
  z <- zs[which.min(ps)]


  out.final <- list(s=s,
                    cp=cp,
                    z=z,
                    mag=mag,
                    p=p.val)

  class(out.final) <- c("list","lacpd.mw")
  attr(out.final,"zs") <- zs
  attr(out.final,"mags") <- mags
  attr(out.final,"ps") <- ps


  return(out.final)
}

print.lacpd.mw <- function(x,round=5,...){
  cat("LACPD Mann-Whitney \n");
  cat("adaptive window:", " ", paste0(x$window),"\n");
  cat("   change-point:", " ", paste0(x$cp),"\n");
  cat("     statistics:",paste0(round(x$z,round)),"\n");
  cat("      magnitude:",paste0(round(x$mag,round)),"\n");
  cat("        p.value:",paste0(round(x$p,round)),"\n");
}
