##' @title Plot of the deviance to derive the confidence interval of of the scale parameter, phi
##' @description This function computes the confidence interval of phi
##' @param obj object of class "SDALGCP" from the the call to function \link{SDALGCPMCML}
##' @param coverage the coverage probability, default is 0.95
##' @param plot logical, to plot the deviance curve. default is TRUE
##' @details This function computes the confidence interval of phi
##' @return return the deviance plot and the corresponding confidence interval of the scale parameter phi
##' @seealso \link{SDALGCPMCML}
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' @importFrom stats loess predict splinefun optimize qchisq 
##' @importFrom graphics abline plot

phiCI <- function(obj, coverage=0.95, plot=TRUE){
  dummy_phiCI <- function(phi, llik.val, coverage=0.95){
    jjj <- cbind(phi, llik.val)
    lo <- stats::loess(llik.val~phi)
    if (length(phi) < 1000){
      x1 <- seq(min(phi),max(phi), length.out = 1000)
    }
    else
    {
      x1 <- seq(min(phi),max(phi), length.out = 2*length(phi))
    }
    y1 <- stats::predict(lo,x1)
    ###############################
    f <- stats::splinefun(x1, y1)
    est.phi <- stats::optimize(function(x) -f(x), lower = min(x1), 
                               upper = max(x1))
    max.log.lik <- -est.phi$objective
    D <- -2*(y1-max.log.lik)
    
    #df <- data.frame(phi=x1, Dev=D)
    #df[which(df$Dev == cut.off), ]
    a <- sign(D - 4)
    b <- which(diff(a) != 0)
    c <- x1[b]
    # d <- c
    # max.phi.est <- jjj[which.max(jjj[,2]),][1]
    # if(is.na(d[1])) d[1] <- min(phi)
    # if(is.na(d[2])) d[2] <- max(phi)
    # 
    # if(d[1] > max.phi.est & d[2] > max.phi.est ) {
    #   c[2] <- min(d)
    #   c[1] <- min(phi)
    # }
    # if(d[1] < max.phi.est & d[2] < max.phi.est ) {
    #   c[1] <- max(d)
    #   c[2] <- max(phi)
    # }
    if(plot){
      graphics::plot(x1, D, xlab=expression(phi), ylab='Deviance', type='l')
      cut.off <- stats::qchisq(coverage, df = 1)
      graphics::abline(h=cut.off, col=2)
      graphics::abline(v=c[1], lty=2, col='blue')
      graphics::abline(v=c[2], lty=2, col='blue')
    }
    return(c(c[1], c[2]))
  }
  out <- dummy_phiCI(obj$all_para[,1], obj$all_para[,2])
  return(out)
}


