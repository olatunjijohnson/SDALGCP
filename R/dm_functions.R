##' @importFrom graphics  abline plot
##' @importFrom methods is
##' @importFrom stats coef dist glm loess median model.frame model.matrix model.offset
##' model.response optimize pnorm predict printCoefmat qchisq qnorm rnorm runif sd splinefun
#if(getRversion() >= "2.15.1") utils::globalVariables(c("j"))
##' @title Generating points inside each polygon in the entire domain using Simple Sequential Inhibition (SSI) process
##' @description This function generate a random point pattern using Simple Sequential Inhibition (SSI) process.
##' @param poly polygon in which to generate the points.
##' @param delta distance between points.
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param pop_shp Optional, The raster of population density map for population weighted approach.
##' @param lambdamax the maximum value of the population density in the polygon.
##' @param pop the population density.
##' @param n optional; the number of points to create in the polygon, if not supplied, it is computed as \eqn{n = rho*|A|*4/(\pi*delta^2)}
##' @param rho Optional, The packing density, default set to 0.55.
##' @param giveup Number of rejected proposals after which the algorithm should terminate.
##' @details This algorithm generates points inside the polygon using Simple Sequential Inhibition (SSI) process.
##' @return It returns a list of the coordinates of the points created in each polygon.
##' @examples 
##' \dontrun{
##' data(PBCshp)
##' My_SSIP(poly=PBCshp@polygons[[1]]@Polygons[[1]]@coords, delta=100)
##' }
##' @importFrom raster extract aggregate
##' @importFrom graphics axis
##' @importFrom spatstat as.owin ppp
##' @importFrom sp SpatialPolygons Polygons Polygon spsample
##' @importFrom splancs areapl csr
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' 
My_SSIP <- function(poly, delta, weighted=FALSE, pop_shp=NULL, lambdamax=NULL, pop=NULL, n= NULL, rho=NULL, giveup=NULL){
  if (weighted==TRUE){
    if (is.null(rho)) rho <- 0.55
    if (is.null(giveup))  giveup <- 1000
    if (is.null(n)) n <- round((rho*splancs::areapl(poly)*4)/ (pi*delta^2))
    if (is.null(pop)) stop("please supply the population total for the region(s) 
                           or if you don't have the population density map just set pop_shp=NULL")
    if (is.null(lambdamax))  stop("please supply the max population for the region(s) 
                                  or if you don't have the population density map just set pop_shp=NULL")
    out.wei <- c()
    fg <- 0
    repeat {
      fg <- fg + 1
      xy<-matrix(splancs::csr(poly,1),1,2)
      vfv <- raster::extract(pop_shp, t(c(xy[1], xy[2])))
      prob2 <- vfv/lambdamax
      wei <- vfv/pop
      if ( fg > 20) {
        wei <- 1/n
        prob2 <- vfv/max(as.numeric(pop_shp@data@values), na.rm = TRUE)
      } 
      if (!is.na(prob2) & prob2 < runif(1)) break
    }
    out.wei <- c(out.wei, wei)
    delsq <- delta*delta
    if ((n * pi * delsq/4 > splancs::areapl(poly))) 
      warning(paste("Window is too small to fit", n, "points", 
                    "at minimum distance", delta))
    len.xy <- 0
    while (len.xy<n) {
      dsq<-0
      k <- 0
      prob <- 0
      while (dsq<delsq | prob < runif(1)) {
        k <- k+1
        ############
        repeat{
          xy.try <- c(splancs::csr(poly,1))
          dsq <- min((xy[,1]-xy.try[1])^2+(xy[,2]-xy.try[2])^2)
          #########rejection sampling
          vdv <- raster::extract(pop_shp, t(c(xy.try[1], xy.try[2])))
          prob <- vdv/lambdamax
          wei <- vdv/pop
          if (!is.na(prob)) break
        }
        #this is to reduce the distance for areas of increased probability of acceptance
        #delsq <- (delta*delta)/exp(prob)
        if (prob<1){
          delsq <- (delta*delta)*(1-prob)
        }else{
          delsq <- (delta*delta)
        }
        ################
        #this to further reduce the computational burden if it is not accepting
        if (k > 20){
          wei <- 1/n
          prob <- vdv/max(as.numeric(pop_shp@data@values), na.rm = TRUE)
        } 
        
        #########termination criteria
        if (k >= giveup) {
          len.xy <- n
          warning(paste("Gave up after", k, "attempts with only",
                        dim(xy)[1], "points placed out of", n))
        }
      }
      if (len.xy ==n) break
      xy<-rbind(xy,xy.try)
      out.wei <- c(out.wei, wei)
      len.xy <- dim(xy)[1]
    }
    rownames(xy) <- 1:nrow(xy)
    colnames(xy) <- c('x', 'y')
    return(list(xy=xy, weight=out.wei/sum(out.wei)))
  }else { #start when the population density is not provided
    if (is.null(rho))  rho <- 0.55
    if (is.null(giveup))  giveup <- 1000
    if (is.null(n)) n <- round((rho*splancs::areapl(poly)*4)/ (pi*delta^2))
    xy <- matrix(splancs::csr(poly,1),1,2)
    delsq <- delta*delta
    if ((n * pi * delsq/4 > splancs::areapl(poly))) 
      warning(paste("Window is too small to fit", n, "points", 
                    "at minimum distance", delta))
    len.xy <- 0
    while (len.xy<n) {
      dsq<-0
      k <- 0
      while (dsq<delsq) {
        k <- k+1
        xy.try <- c(splancs::csr(poly,1))
        dsq <- min((xy[,1]-xy.try[1])^2+(xy[,2]-xy.try[2])^2)
        #########termination criteria
        if (k >= giveup) {
          len.xy <- n
          warning(paste("Gave up after", k, "attempts with only",
                        dim(xy)[1], "points placed out of", n))
        }
      }
      if (len.xy ==n) break
      xy<-rbind(xy,xy.try)
      len.xy <- dim(xy)[1]
    }
    rownames(xy) <- 1:nrow(xy)
    colnames(xy) <- c('x', 'y')
    return(list(xy=xy))
  }
}
####################################
##' @title Generating points inside each polygon in the entire domain using a uniform sampling or completely spatial random sampling.
##' @description This function generate a random point pattern using a uniform sampling or completely spatial random sampling.
##' @param poly polygon in which to generate the points
##' @param delta distance between points
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param pop_shp Optional, The raster of population density map for population weighted approach
##' @param lambdamax the maximum value of the population density in the polygon
##' @param pop the population density.
##' @param n optional; the number of points to create in the polygon, if not supplied, it is computed as \eqn{n = rho*|A|*4/(\pi*delta^2)}
##' @param rho Optional, The packing density, default set to 0.55
##' @param giveup Number of rejected proposals after which the algorithm should terminate.
##' @param bound Spatial object; the boundary of the polygon
##' @details This algorithm generates points inside the polygon using a uniform sampling or completely spatial random sampling.
##' @return It returns a list of the coordinates of the points created in each polygon.
##' @examples 
##' \dontrun{
##' data(PBCshp)
##' My_UNIF(poly=PBCshp@polygons[[1]]@Polygons[[1]]@coords, delta=100)
##' }
##' @importFrom raster extract aggregate
##' @importFrom graphics axis
##' @importFrom spatstat as.owin ppp
##' @importFrom sp SpatialPolygons Polygons Polygon spsample
##' @importFrom splancs areapl csr
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' 
My_UNIF <- function(poly, delta, weighted=FALSE, pop_shp=NULL, lambdamax=NULL, pop=NULL, n= NULL, rho=NULL, giveup=NULL, 
                    bound=NULL){
  if (weighted==TRUE){
    if (is.null(rho)) rho <- 0.55
    if (is.null(giveup))  giveup <- 1000
    if (is.null(n)) n <- round((rho*splancs::areapl(poly)*4)/ (pi*delta^2))
    if (is.null(pop)) stop("please supply the population total for the region(s) 
                           or if you don't have the population density map just set pop_shp=NULL")
    if (is.null(lambdamax))  stop("please supply the max population for the region(s) 
                                  or if you don't have the population density map just set pop_shp=NULL")
    out.wei <- c()
    fg <- 0
    repeat {
      fg <- fg + 1
      xy<-matrix(splancs::csr(poly,1),1,2)
      vfv <- raster::extract(pop_shp, t(c(xy[1], xy[2])))
      prob2 <- vfv/lambdamax
      wei <- vfv/pop
      if ( fg > 20) {
        wei <- 1/n
        prob2 <- vfv/max(as.numeric(pop_shp@data@values), na.rm = TRUE)
      } 
      if (!is.na(prob2) & prob2 < runif(1)) break
    }
    out.wei <- c(out.wei, wei)
    delsq <- delta*delta
    if ((n * pi * delsq/4 > splancs::areapl(poly))) 
      warning(paste("Window is too small to fit", n, "points", 
                    "at minimum distance", delta))
    len.xy <- 0
    while (len.xy<n) {
      k <- 0
      prob <- 0
      while (prob < runif(1)) {
        k <- k+1
        ############
        repeat{
          xy.try <- c(splancs::csr(poly,1))
          #########rejection sampling
          vdv <- raster::extract(pop_shp, t(c(xy.try[1], xy.try[2])))
          prob <- vdv/lambdamax
          wei <- vdv/pop
          if (!is.na(prob)) break
        }
        if (k > 20){
          wei <- 1/n
          prob <- vdv/max(as.numeric(pop_shp@data@values), na.rm = TRUE)
        } 
        #########termination criteria
        if (k >= giveup) {
          len.xy <- n
          warning(paste("Gave up after", k, "attempts with only",
                        dim(xy)[1], "points placed out of", n))
        }
      }
      if (len.xy ==n) break
      xy<-rbind(xy,xy.try)
      out.wei <- c(out.wei, wei)
      len.xy <- dim(xy)[1]
    }
    rownames(xy) <- 1:nrow(xy)
    colnames(xy) <- c('x', 'y')
    return(list(xy=xy, weight=out.wei/sum(out.wei)))
  }else { #start when the population density is not provided
    if (is.null(rho)) rho <- 0.55
    if (is.null(giveup)) giveup <- 1000
    if (is.null(n)) n <- round((rho*splancs::areapl(poly)*4)/ (pi*delta^2))
    xy<-matrix(splancs::csr(poly,1),1,2)
    delsq <- delta*delta
    if ((n * pi * delsq/4 > splancs::areapl(poly))) 
      warning(paste("Window is too small to fit", n, "points", 
                    "at minimum distance", delta))
    xy <- sp::spsample(bound, cellsize=delta, type = "random")
    rownames(xy) <- 1:nrow(xy)
    colnames(xy) <- c('x', 'y')
    return(list(xy=xy))
  }
}

###########################
##' @title Generating points inside each polygon in the entire domain using a regular (systematically aligned) sampling.
##' @description This function generate a random point pattern using a regular (systematically aligned) sampling.
##' @param poly polygon in which to generate the points
##' @param delta distance between points
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param pop_shp Optional, The raster of population density map for population weighted approach.
##' @param lambdamax the maximum value of the population density in the polygon
##' @param pop the population density.
##' @param n optional; the number of points to create in the polygon, if not supplied, it is computed as \eqn{n = rho*|A|*4/(\pi*delta^2)}
##' @param rho Optional, The packing density, default set to 0.55
##' @param giveup Number of rejected proposals after which the algorithm should terminate.
##' @param bound Spatial object; the boundary of the polygon
##' @details This algorithm generates points inside the polygon using a regular (systematically aligned) sampling.
##' @return It returns a list of the coordinates of the points created in each polygon.
##' @examples 
##' \dontrun{
##' data(PBCshp)
##' My_REG(poly=PBCshp@polygons[[1]]@Polygons[[1]]@coords, delta=100)
##' }
##' @importFrom raster extract aggregate
##' @importFrom graphics axis
##' @importFrom spatstat as.owin ppp
##' @importFrom sp SpatialPolygons Polygons Polygon spsample
##' @importFrom splancs areapl csr
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' 
My_REG <- function(poly, delta, weighted=FALSE, pop_shp=NULL, lambdamax=NULL, pop=NULL, n= NULL, rho=NULL, giveup=NULL,
                   bound=NULL){
  if (weighted==TRUE){
    stop("There is no option to use regular grid with population density. Please use uniform or SSIP")
  }else { #start when the population density is not provided
    if (is.null(rho)) rho <- 0.55
    if (is.null(giveup))  giveup <- 1000
    if (is.null(n)) n <- round((rho*splancs::areapl(poly)*4)/ (pi*delta^2))
    xy<-matrix(splancs::csr(poly,1),1,2)
    delsq <- delta*delta
    if ((n * pi * delsq/4 > splancs::areapl(poly))) 
      warning(paste("Window is too small to fit", n, "points", 
                    "at minimum distance", delta))
    xy <- sp::spsample(bound, cellsize=delta, type = "regular")
    rownames(xy) <- 1:nrow(xy)
    colnames(xy) <- c('x', 'y')
    return(list(xy=xy))
  }
}
########################################################
##' @title The wrapper function for the three methods to be used to create points inside the polygon
##' @description This function generate a random point pattern using Simple Sequential Inhibition (SSI) process, uniform sampling and regular grid point.
##' @param my_shp A SpatialPolygons orSpatialPolygonsDataFrame  object containing the polygons (i.e each regions).
##' @param delta distance between points
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param lambdamax the maximum value of the population density in the polygon.
##' @param pop the population density.
##' @param pop_shp Optional, The raster of population density map for population weighted approach.
##' @param n optional; the number of points to create in the polygon, if not supplied, it is computed as \eqn{n = rho*|A|*4/(\pi*delta^2)}
##' @param method To specify which method to use to sample the points, the options are 1 for Simple Sequential Inhibition (SSI) process, 2 for Uniform sampling and 3 for regular grid. 1 is the default
##' @param plot To display the plot of the points inside the polygon, default to TRUE.
##' @param rho Optional, The packing density, default set to 0.55
##' @param giveup Number of rejected proposals after which the algorithm should terminate.
##' @details This algorithm generates points inside the polygon using three algorithms specified in the method. 
##' @return It returns a list of the coordinates of the points created in each polygon and it has an associated attribute weighted which is either TRUE or FALSE to indicate if the population density is used or not.
##' @examples 
##' \dontrun{
##' data(PBCshp)
##' create_points(my_shp=PBCshp, delta=100)
##' }
##' @importFrom raster extract aggregate
##' @importFrom graphics axis
##' @importFrom spatstat as.owin ppp
##' @importFrom sp SpatialPolygons Polygons Polygon spsample
##' @importFrom splancs areapl csr
##' @importFrom maptools as.owin.SpatialPolygons
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
create_points <- function(my_shp, delta, weighted=FALSE, lambdamax=NULL, pop=NULL, pop_shp=NULL, n=NULL, method=1,
                          plot=FALSE, rho=NULL, giveup=NULL){
  #################################################
  if (method==1){
    SSI=TRUE
  }else if (method==2){
    uniform=TRUE
  }else{
    regular=TRUE
  }
  ###############################################
  if (class(my_shp)[1]=="Polygon"){
    bound <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(xrange=my_shp@coords[,1], 
                                                                          yrange=my_shp@coords[,2]))), ID = "a")))
    poly <- bound@polygons[[1]]@Polygons[[1]]@coords
    #r.shp <- bbox(my_shp)
  } else{
    bound <- raster::aggregate(my_shp)
    poly <- bound@polygons[[1]]@Polygons[[2]]@coords
    #r.shp <- bbox(bound)
  }
  #check lambdamax if it is supposed to be set at this value or to stop
  #if (is.null(lambdamax)) lambdamax <- max(as.numeric(pop_shp@data@values), na.rm = TRUE)
  if (SSI==TRUE){
    xycand <- My_SSIP(poly=poly, delta=delta, pop_shp=pop_shp, lambdamax=lambdamax, pop=pop,
                      n = n, rho=rho, giveup=giveup, weighted=weighted)
    
  } else if (uniform==TRUE){
    xycand <- My_UNIF(poly=poly, delta=delta, pop_shp=pop_shp, lambdamax=lambdamax, pop=pop,
                      n = n, rho=rho, giveup=giveup, weighted=weighted, bound=bound)
  } else{
    xycand <- My_REG(poly=poly, delta=delta, pop_shp=pop_shp, lambdamax=lambdamax, pop=pop,
                     n = n, rho=rho, giveup=giveup, weighted=weighted, bound=bound)
  }
  
  if (plot==TRUE) {
    if (class(my_shp)[1]=="Polygon"){
      cases.sim <- spatstat::ppp(x = xycand$xy[,1], y = xycand$xy[,2], window = maptools::as.owin.SpatialPolygons(bound))
    }else{
      cases.sim <- spatstat::ppp(x = xycand$xy[,1], y = xycand$xy[,2], window = spatstat::as.owin(bound))
    }
    plot(cases.sim)
    graphics::axis(1)
    graphics::axis(2)
  }
  return(xycand)
}

#########################################
##' @title Generating points inside each polygon in the entire domain
##' @description This function generate a random point pattern using Simple Sequential Inhibition (SSI) process, uniform sampling and regular grid point.
##' @param my_shp A SpatialPolygons orSpatialPolygonsDataFrame  object containing the polygons (i.e each regions).
##' @param delta distance between points
##' @param pop_shp Optional, The raster of population density map for population weighted approach
##' @param rho Optional, The packing density, default set to 0.55
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param plot To display the plot of the points inside the polygon, default to TRUE
##' @param method To specify which method to use to sample the points, the options are 1 for Simple Sequential Inhibition (SSI) process, 2 for Uniform sampling and 3 for regular grid. 1 is the default
##' @details This algorithm generates points inside the polygon using three algorithms specified in the method. 
##' @return It returns a list of the coordinates of the points created in each polygon and it has an associated attribute weighted which is either TRUE or FALSE to indicate if the population density is used or not.
##' @importFrom raster extract aggregate
##' @importFrom graphics axis
##' @importFrom spatstat as.owin ppp
##' @importFrom sp SpatialPolygons Polygons Polygon spsample
##' @importFrom splancs areapl csr
##' @examples 
##' \dontrun{
##' data(PBCshp)
##' SDALGCPpolygonpoints(my_shp=PBCshp, delta=100)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
SDALGCPpolygonpoints <- function(my_shp, delta, method=1, pop_shp=NULL,  weighted=FALSE, rho=NULL, plot=TRUE){
  if(weighted == TRUE & is.null(pop_shp)) stop('please insert the raster file of the population density or change argument weights==FALSE if you do not plan to use population density')
  if(weighted==FALSE){
    my_list <- list()
    for (i in 1:length(my_shp)){
      my_list[[i]] <- create_points(my_shp = my_shp@polygons[[i]]@Polygons[[1]], 
                                    pop_shp = pop_shp, delta=delta, method=method,
                                    plot=plot, lambdamax=NULL, pop=NULL, rho=NULL, weighted=weighted)
      cat('creating points inside region', i, 'out of', length(my_shp), 'regions', '\n')
    }
    attr(my_list, 'weighted') <- FALSE
    attr(my_list, 'my_shp') <- my_shp
    return(my_list)
  } else{
    cat("Extracting the population density for each polygon")
    pop_lsoa <- raster::extract(pop_shp, my_shp, weights=weighted,normalizeWeights=F)
    summ.mat <- function(my_mat) {
      #my_answer <- my_mat[,1] %*% my_mat[,2]
      my_answer <- sum(my_mat[,1]*my_mat[,2])
      return(round(my_answer))
    }
    my_pop_lsoa <- unlist(lapply(pop_lsoa, FUN = summ.mat))
    max.mat <- function(my_mat){
      my_answer <- max(my_mat[,1])
      return(round(my_answer))
    }
    my_pop_lsoa_max <- unlist(lapply(pop_lsoa,FUN = max.mat))
    my_list <- list()
    for (i in 1:length(my_shp)){
      my_list[[i]] <- create_points(my_shp = my_shp@polygons[[i]]@Polygons[[1]], 
                                    pop_shp = pop_shp, delta=delta, method=method, rho=rho,
                                    plot=plot, lambdamax=my_pop_lsoa_max[i], pop=my_pop_lsoa[i], weighted=weighted)
      cat('creating points inside region', i, 'out of', length(my_shp), 'regions', '\n')
    }
    attr(my_list, 'weighted') <- TRUE
    attr(my_list, 'my_shp') <- my_shp
    return(my_list)
  }
}

###############################################################
##' @title Pre-compute the covariance matrix
##' @description This function precomputes the covariance matrix 
##' @param S.coord The list of the coordinates of the points created in each polygon
##' @param phi The vector of the scale parameter
##' @details This function precompute the covariance matrix using the exponential covariance function
##' @return return an array of the covariance matrix for the specified vectors of phi
##' @examples
##' \dontrun{
##' data(PBCshp)
##' S.coord <- SDALGCPpolygonpoints(my_shp=PBCshp, delta=100)
##' phi <- seq(1, 10, 1)
##' pre.matrix <- precomputeCorrMatrix(S.coord=S.coord, phi=phi)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom utils setTxtProgressBar txtProgressBar
##' @seealso \code{\link{SDALGCPpolygonpoints}}, 
##' @export

precomputeCorrMatrix <- function(S.coord, phi){
  weight=attr(S.coord, 'weighted')
  n.distr <- length(S.coord)
  n.phi <- length(phi)
  cat("\n Start precomputing the correlation matrix!\n")
  pb = utils::txtProgressBar(min = 0, max = n.distr, initial = 0) 
  R= array(NA, dim = c(n.distr, n.distr, n.phi))
  if (weight==TRUE){
    for (i in 1:n.distr){
      utils::setTxtProgressBar(pb,i, label=paste( round(i/n.distr*100, 0), "% done"))
      Sys.sleep(0.01)
      for (j in i:n.distr){
        U <- as.matrix(pdist::pdist(as.matrix(S.coord[[i]]$xy), as.matrix(S.coord[[j]]$xy)))
        W <- outer(S.coord[[i]]$weight, S.coord[[j]]$weight, '*')
        for (k in 1:n.phi){
          R[i,j, k] <- R[j,i, k] <- sum(W*exp(-U/phi[k]))
        }
      }
    }
    attr(R, 'weighted') <- TRUE
    attr(R, 'my_shp') <-   attr(S.coord, 'my_shp')
    attr(R, 'S_coord') <-   S.coord   
  }else{
    for (i in 1:n.distr){
      utils::setTxtProgressBar(pb,i, label=paste( round(i/n.distr*100, 0), "% done"))
      Sys.sleep(0.01)
      for (j in i:n.distr){
        U <- as.matrix(pdist::pdist(as.matrix(S.coord[[i]]$xy), as.matrix(S.coord[[j]]$xy)))
        for (k in 1:n.phi){
          R[i,j, k] <- R[j,i, k] <- mean(exp(-U/phi[k]))
        }
      }
    }
    attr(R, 'weighted') <- FALSE
    attr(R, 'my_shp') <-   attr(S.coord, 'my_shp')
    attr(R, 'S_coord') <-   S.coord   
  }
  close(pb)
  cat("\n Done precomputing the correlation matrix!\n")
  return(list(R=R, phi=phi))
}

##################################################################################
##' @title Monte Carlo Maximum Likelihood estimation for the Poisson model of aggregated data
##' @description This function performs Monte Carlo maximum likelihood (MCML) estimation for the geostatistical Poisson model with log link function.
##' @param y the data
##' @param D the design matrix
##' @param m the offset term
##' @param corr the correlation matrix from exponential correlation function
##' @param par0 the initial parameter of the fixed effects beta, the variance sigmasq and the scale parameter phi, specified as c(beta, sigma2, phi)
##' @param control.mcmc output from \code{\link{controlmcmcSDA}}.
##' @param S.sim the posterior sample of the linear predictor given the initial parameters
##' @param Denominator the value of the denominator of the likelihood
##' @details The function helps to obtain the MCML estimate for a given value of correlation matrix, i.e for a given value of the scale parameter phi.
##' @return \code{estimate}: estimates of the model parameters; beta's and with sigma2 on the log scale
##' @return \code{covariance}: covariance matrix of the MCML estimates.
##' @return \code{log.lik}: maximum value of the log-likelihood.
##' @return \code{S}: the linear predictor given the initial parameter
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @references Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. Journal of Statistical Software, 78(8), 1-29. doi:10.18637/jss.v078.i08.
##' @references Christensen, O. F. (2004). Monte carlo maximum likelihood in model-based geostatistics. Journal of Computational and Graphical Statistics 13, 702-718.
##' @importFrom stats nlminb 
##' @importFrom PrevMap Laplace.sampling
##' @keywords internal
##' @seealso \code{\link{controlmcmcSDA}}

Aggregated_poisson_log_MCML <- function(y, D, m, corr, par0, control.mcmc, S.sim,
                                        Denominator) {
  n <- length(y)
  p <- ncol(D)
  
  
  R.inv <- solve(corr)
  ldetR <- determinant(corr)$modulus
  
  #likelihood
  Log.Joint.dens.S.Y <- function(S,val) {
    llik <- sum(y*S-m*exp(S))
    diff.S <- S-val$mu
    AAA <-    t(diff.S)%*%R.inv%*%(diff.S)
    return(-0.5*(n*log(val$sigma2)+ ldetR+
                   AAA/val$sigma2)+ llik)
  }
  
  
  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    val$mu <- as.numeric(D%*%beta)
    val$sigma2 <- exp(par[p+1])
    
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,],val)))
  }
  
  Monte.Carlo.Log.Lik <- function(par) {
    log(mean(exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)))  
  }
  
  
  grad.Monte.Carlo.Log.Lik <- function(par){
    beta <- par[1:p]
    D.beta <- D%*%beta
    sigma2 <- exp(par[p+1])
    
    First.deriv.S.param <- function(S){
      
      diff.S <- S-D.beta
      AAA <- t(diff.S)%*%R.inv%*%diff.S
      grad.beta <-  t(D)%*%R.inv%*%(diff.S)/sigma2
      
      grad.log.sigma2 <- (-n/(2*sigma2)+0.5*AAA/(sigma2^2))*sigma2
      
      der.par <- c(grad.beta, grad.log.sigma2)
      return(der.par)
    }
    
    
    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio            
    
    
    cumm <- rep(0,length(par))
    for(i in 1:(dim(S.sim)[1])) {
      full.deriv <- part.deriv[i]*First.deriv.S.param(S.sim[i,])
      cumm <- cumm + full.deriv
    }
    return(cumm)
  }
  
  
  
  #The second derivative of the Monte Carlo approximation
  hess.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    mu <- D%*%beta
    sigma2 <- exp(par[p+1])
    
    H <- matrix(0,nrow=length(par),ncol=length(par))
    H[1:p,1:p] <- (-t(D)%*%R.inv%*%D)/sigma2
    
    Second.deriv.S.param <- function(S, part.deriv) {
      
      diff.S <- S-mu
      q.f <- t(diff.S)%*%R.inv%*%diff.S
      
      grad.beta <-  t(D)%*%R.inv%*%(diff.S)/sigma2
      
      grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
      
      der.par <- c(grad.beta, grad.log.sigma2)
      
      
      H[1:p,p+1] <- H[p+1,1:p] <- -t(D)%*%R.inv%*%(diff.S)/sigma2
      
      
      H[p+1,p+1] <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
        grad.log.sigma2
      
      out <- list()
      out$first.term<- part.deriv*(der.par%*%t(der.par)+H)
      out$grad <- der.par*part.deriv
      out
    }
    
    
    
    ratio <- exp(Num.Monte.Carlo.Log.Lik(par)-Denominator)
    sum.ratio <- sum(ratio)
    part.deriv <- ratio/sum.ratio
    
    
    last.term <- rep(0,length(par))
    cumm <- matrix(0,length(par),length(par))
    for(i in 1:(dim(S.sim)[1])) {
      Hess <- Second.deriv.S.param(S.sim[i,], part.deriv[i])
      last.term <- last.term + Hess$grad
      cumm <- cumm + Hess$first.term
    }
    return(cumm-last.term%*%t(last.term))
  }
  
  new.par <- par0[-length(par0)]
  new.par[(p+1)] <- log(new.par[(p+1)])
  
  output <- list()
  
  result <- stats::nlminb(new.par,function(x) -Monte.Carlo.Log.Lik(x),
                          function(x) -grad.Monte.Carlo.Log.Lik(x),
                          function(x) -hess.Monte.Carlo.Log.Lik(x),control=list(trace=0))
  #i can change trace =0, so that it doesn't print result
  output$estimate <- result$par
  output$covariance <- solve(-hess.Monte.Carlo.Log.Lik(result$par))
  output$value <- -result$objective
  output$S <- S.sim
  names(output$estimate)[1:p] <- colnames(D)
  names(output$estimate)[(p+1)] <- c("sigma^2") #note that it stored log(sigma^2)
  rownames(output$covariance) <- colnames(output$covariance) <- names(output$estimate)
  return(output)
}

######################################################################
##' @title Parameter estimation for SDA-LGCP 
##' @description This function provides the maximum likelihood estimation of the parameter given the precomputed correlation matrices for different values of scale parameter, phi
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param data  data frame containing the variables in the model.
##' @param corr the array of the precomputed correlation matrix for each value of the scale parameter.
##' @param par0 the initial parameter of the fixed effects beta, the variance sigmasq and the scale parameter phi, specified as c(beta, sigma2, phi)
##' @param control.mcmc list from PrevMap package to define the burnin, thining, the number of iteration and the turning parameters see \code{\link{controlmcmcSDA}}.
##' @param plot_profile logical; if TRUE the profile-likelihood is plotted. default is FALSE
##' @details This function performs parameter estimation for a SDA-LGCP Model
##' \bold{Monte Carlo Maximum likelihood.}
##' The Monte Carlo maximum likelihood method uses conditional simulation from the distribution of the random effect \eqn{T(x) = d(x)'\beta+S(x)} given the data \code{y}, in order to approximate the high-dimensiional intractable integral given by the likelihood function. The resulting approximation of the likelihood is then maximized by a numerical optimization algorithm which uses analytic epression for computation of the gradient vector and Hessian matrix. The functions used for numerical optimization are \code{\link{nlminb}} 
##' @return An object of class "PrevMap".
##' The function \code{\link{summary.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{D}: matrix of covariates.
##' @return \code{y}: the count, response observations.
##' @return \code{m}: offset
##' @return \code{beta_opt}: estimates of the fixed effects of the model.
##' @return \code{sigma2_opt}: estimates of the variance of the Gaussian process.
##' @return \code{phi_opt}: estimates of the scale parameter phi of the Gaussian process.
##' @return \code{cov}: covariance matrix of the MCML estimates.
##' @return \code{Sigma_mat_opt }: covariance matrix of the Gaussian process that corresponds to the optimal value
##' @return \code{llike_val_opt}: maximum value of the log-likelihood.
##' @return \code{mu}: mean of the linear predictor
##' @return \code{all_para}: the entire estimates for the different values of phi.
##' @return \code{all_cov}: the entire covariance matrix of the estimates for the different values of phi.
##' @return \code{par0}: the initial parameter of the fixed effects beta and the variance sigmasq used in the estimation
##' @return \code{control.mcmc}:  the burnin, thining, the number of iteration and the turning parameters used see \code{\link{controlmcmcSDA}}.
##' @return \code{S}: the linear predictor given the initial parameter
##' @return \code{call}: the matched call.
##' @examples
##' \dontrun{
##' data(PBCshp)
##' S.coord <- SDALGCPpolygonpoints(my_shp=PBCshp, delta=100)
##' phi <- seq(1, 10, 1)
##' pre.matrix <- precomputeCorrMatrix(S.coord=S.coord, phi=phi)
##' data <- PBCshp@data
##' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
##' Environment +  offset(log(pop))
##' SDAParaEst(formula=FORM , data=data, corr=pre.matrix)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom pbapply pblapply
##' @references Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. Journal of Statistical Software, 78(8), 1-29. doi:10.18637/jss.v078.i08.
##' @references Christensen, O. F. (2004). Monte carlo maximum likelihood in model-based geostatistics. Journal of Computational and Graphical Statistics 13, 702-718.
##' @seealso \link{Aggregated_poisson_log_MCML}, \code{\link{Laplace.sampling}}
##' @export

SDAParaEst <- function(formula, data, corr, par0=NULL, control.mcmc=NULL, plot_profile=FALSE){
  cat("\n Now preparing for parameter estimation!\n")
  mf <- model.frame(formula=formula,data=data)
  y <- as.numeric(model.response(mf))
  D <- model.matrix(attr(mf,"terms"), data=data)
  n <- length(y)
  p <- ncol(D)
  if(any(startsWith(names(mf), 'offset')==TRUE)) {
    m <-  exp(model.offset(mf))
  } else {
    m <- rep(1, n)
  }
  if(is.null(par0)) {
    phi <- as.numeric(corr$phi)                
    n.phi <- length(phi)
    R <- corr$R
    model <- glm(formula, family="poisson", data=data)
    beta.start <-coef(model)
    sigma2.start <- mean(model$residuals^2)
    phi.start <- median(phi)
    par0 <- c(beta.start, sigma2.start, phi.start)
    whichmedian <- function(x) which.min(abs(x - median(x)))
    corr0 <- R[,,whichmedian(phi)]
  }else{
    phi <- as.numeric(corr$phi)
    phi <- phi[-(length(phi))]
    n.phi <- length(phi)
    R <- corr$R[,,(-(n.phi+1))]
    corr0 <- R[,,(n.phi+1)]
  }
  if(any(par0[-(1:p)] <= 0)) stop("the covariance parameters in 'par0' must be positive.")
  if(is.null(control.mcmc)) control.mcmc <- list(n.sim = 10000, burnin = 2000, thin= 8, h=1.65/(n^(1/6)),
                                                 c1.h = 0.01, c2.h = 1e-04)

  #######################################MCMC
  #initial values
  beta0 <- par0[1:p]
  mu0 <- as.numeric(D%*%beta0)
  sigma2.0 <- par0[p+1]
  Sigma0 <- sigma2.0 * corr0
  cat("\n Simulating the linear predictor given the initial parameter \n")
  repeat {
    S.sim.res <- tryCatch(PrevMap::Laplace.sampling(mu=mu0, Sigma=Sigma0, y=y, units.m=m, 
                                                    control.mcmc=control.mcmc,
                                                    plot.correlogram=FALSE,messages=FALSE,
                                                    poisson.llik=TRUE), error=identity)
    if (!is(S.sim.res, "error"))
      break
  }
  S.sim <- S.sim.res$samples
  R.inv0 <- solve(corr0)
  ldetR0 <- determinant(corr0)$modulus
  ################# compute the denominator
  Log.Joint.dens.S.Y <- function(S,val) {
    llik <- sum(y*S-m*exp(S))
    diff.S <- S-val$mu
    AAA <-    t(diff.S)%*%R.inv0%*%(diff.S)
    return(-0.5*(n*log(val$sigma2)+ ldetR0+
                   AAA/val$sigma2)+ llik)
  }
  
  #it computes the density of S for each sample of S
  Num.Monte.Carlo.Log.Lik <- function(par) {
    beta <- par[1:p]
    val <- list()
    val$mu <- as.numeric(D%*%beta)
    val$sigma2 <- exp(par[p+1])
    
    return(sapply(1:(dim(S.sim)[1]),function(i) Log.Joint.dens.S.Y(S.sim[i,],val)))
  }
  Den.Monte.Carlo.Log.Lik <- Num.Monte.Carlo.Log.Lik(c(beta0,log(sigma2.0)))
  
  
  ######################################
  func <- function(x){
    result <- Aggregated_poisson_log_MCML(y=y, D=D, m=m, corr= R[,,x], par0=par0, 
                                          control.mcmc=control.mcmc, S.sim=S.sim, 
                                          Denominator = Den.Monte.Carlo.Log.Lik)
    result$estimate[p+1] <- exp(result$estimate[p+1])
    return(list(par=c(phi[x], result$value, as.numeric(result$estimate)), cov=result$covariance))
  }
  cat("\n Now estimating the parameter \n")
  ress <- pbapply::pblapply(1:n.phi, FUN = func)
  output <- as.data.frame(do.call('rbind', lapply(ress, function(x) x$par)))
  output2 <-  lapply(ress, function(x) x$cov)
  ########to get predictors names
  mt <- attr(mf, "terms")
  predictorsnames <- c("(intercept)", attr(mt, "term.labels"))
  ##########
  colnames(output) <- c('phi', 'value', predictorsnames, 'sigma2')
  #i need to redo the col name when par0 is specified 
  if (plot_profile==TRUE) plot(output[,1], output[,2], type='l', ylab='loglik', xlab='phi')
  max.ind <- which.max(output[,'value'])
  max.res=output[max.ind,]
  colnames(max.res) <- c('phi', 'value', predictorsnames, 'sigma2')
  cov.max.res <- output2[[max.ind]]
  out <- list()
  out$D <- D
  out$y <- y
  out$m <- m
  out$beta_opt <- as.numeric(max.res[predictorsnames])
  out$sigma2_opt <- as.numeric(max.res['sigma2'])
  out$phi_opt <- as.numeric(max.res['phi'])
  out$cov <- cov.max.res
  out$Sigma_mat_opt <- out$sigma2_opt*R[,,which.max(output[,'value'])]
  out$llike_val_opt <- as.numeric(max.res['value'])
  out$mu <- D%*%out$beta_opt
  out$all_para <- output
  out$all_cov <- output2
  out$par0 <- par0
  out$control.mcmc <- control.mcmc
  out$S <- S.sim
  out$call <- match.call()
  attr(out, 'weighted') <- attr(corr$R, 'weighted')
  attr(out, 'my_shp') <- attr(corr$R, 'my_shp')
  attr(out, 'S_coord') <- attr(corr$R, 'S_coord')
  class(out) <- "SDALGCP"
  return(out)
}
######################
#######################################################################################
##' @title Spatial Discrete predictions of the relative risk using plug-in of MCML estimates
##' @description This function performs spatial discrete prediction for a fixed the model parameters at the Monte Carlo maximum likelihood estimates of a SDA-LGCP model.
##' @param para_est an object of class "SDALGCP" obtained as a result of a call to \code{\link{SDAParaEst}}.
##' @param control.mcmc output from \code{\link{controlmcmcSDA}}, if not provided, it uses the values used for the parameter estimation
##' @param divisor optional, if the coordinate of the shapefile is rescaled, default is 1
##' @param plot.correlogram logical; if plot.correlogram=TRUE the autocorrelation plot of the conditional simulations is displayed.
##' @param messages logical; if messages=TRUE then status messages are printed on the screen (or output device) while the function is running. Default is messages=TRUE.
##' @details The function returns the region-specific incidence and the covariate adjusted relative risk exp(S(A))
##' @return It returns the following list
##' @return S.draw: the samples of the linear predictor \eqn{d(x)'\beta + S(A)}
##' @return incidence: the region-specific incidence
##' @return SEincidence: Standard errror of the region-specific incidence
##' @return CovRR: the prediction of the relative risk
##' @return SECovRR: the standard error of the relative risk
##' @examples
##' \dontrun{
##' data(PBCshp)
##' S.coord <- SDALGCPpolygonpoints(my_shp=PBCshp, delta=100)
##' phi <- seq(1, 10, 1)
##' pre.matrix <- precomputeCorrMatrix(S.coord=S.coord, phi=phi)
##' data <- PBCshp@data
##' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
##' Environment +  offset(log(pop))
##' est <- SDAParaEst(formula=FORM , data=data, corr=pre.matrix)
##' SDADiscretePred(para_est=est)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
##' @importFrom sp spsample
##' @importFrom Matrix solve chol
SDADiscretePred <- function(para_est, control.mcmc=NULL,
                            divisor=1, plot.correlogram=FALSE, messages=TRUE){
  out <- list()
  my_shp <- attr(para_est, 'my_shp')
  beta <- para_est$beta_opt 
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt*median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  Sigma0 <- para_est$Sigma_mat_opt
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  S.sim.res <- PrevMap::Laplace.sampling(mu=mu0, Sigma=Sigma0, y=y, 
                                         units.m=m, control.mcmc = control.mcmc,
                                         plot.correlogram=plot.correlogram, messages=messages,
                                         poisson.llik=TRUE)
  S.sim <- S.sim.res$samples
  n.sim <- dim(S.sim)[1]
  my_shp$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]-mu0), 1, mean))
  my_shp$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,]-mu0)), 1, sd)
  my_shp$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]), 1, mean))
  my_shp$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,])), 1, sd)
  out$incidence <- my_shp$pMean_RR
  out$SEincidence <- my_shp$pSD_RR
  out$CovRR <- my_shp$pMean_ARR
  out$SECovRR <- my_shp$pSD_ARR
  out$S.draw <- S.sim
  out$my_shp <- my_shp
  out$para_est <- para_est
  out$call <- match.call()
  attr(out, 'weighted') <- attr(para_est, 'weighted')
  class(out) <- "SDALGCP"
  return(out)
}
#################################################
##' @title Spatial continuous predictions of the relative risk using plug-in of MCML estimates
##' @description This function performs spatial continuous prediction, fixing the model parameters at the Monte Carlo maximum likelihood estimates of a SDA-LGCP model.
##' @param para_est an object of class "SDALGCP" obtained as a result of a call to \code{\link{SDAParaEst}}.
##' @param cellsize the size of the computational grid 
##' @param control.mcmc output from \code{\link{controlmcmcSDA}}, if not provided, it uses the values used for the parameter estimation
##' @param pred.loc optional, the dataframe of the predictive grid.
##' @param divisor optional, the value to use to convert the dimension of the polygon, default is 1 which implies no conversion
##' @param plot.correlogram ogical; if plot.correlogram=TRUE the autocorrelation plot of the conditional simulations is displayed.
##' @param messages logical; if messages=TRUE then status messages are printed on the screen (or output device) while the function is running. Default is messages=TRUE.
##' @param parallel to parallelize some part of the function.
##' @details The function returns the prediction of the relative risk exp(S(x))
##' @return pred.draw: the samples of the prediction
##' @return pred: the prediction of the relative risk
##' @return predSD: the standard error of the prediction
##' @return Pred.loc: The cordinates of the predictive locations
##' @examples
##' \dontrun{
##' data(PBCshp)
##' S.coord <- SDALGCPpolygonpoints(my_shp=PBCshp, delta=100)
##' phi <- seq(1, 10, 1)
##' pre.matrix <- precomputeCorrMatrix(S.coord=S.coord, phi=phi)
##' data <- PBCshp@data
##' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
##' Environment +  offset(log(pop))
##' est <- SDAParaEst(formula=FORM , data=data, corr=pre.matrix)
##' SDAContinuousPred(para_est=est, cellsize=300)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom sp spsample coordinates
##' @importFrom Matrix solve chol
##' @importFrom pdist pdist 
##' @importFrom stats median
##' @export

SDAContinuousPred <- function(para_est, cellsize, control.mcmc=NULL, pred.loc=NULL,
                              divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE){
  out <- list()
  my_shp <- attr(para_est, 'my_shp')
  weight <- attr(para_est, 'weighted')
  S.coord <- attr(para_est, 'S_coord')
  beta <- para_est$beta_opt 
  mu0 <- para_est$mu
  sigma2 <- para_est$sigma2_opt*median(diag(para_est$Sigma_mat_opt))
  phi <-  para_est$phi_opt
  Sigma0 <- para_est$Sigma_mat_opt
  if (is.null(control.mcmc)) control.mcmc <- para_est$control.mcmc
  m <- para_est$m
  y <- para_est$y
  S.sim.res <- PrevMap::Laplace.sampling(mu=mu0, Sigma=Sigma0, y=y, 
                                         units.m=m, control.mcmc = control.mcmc,
                                         plot.correlogram=plot.correlogram, messages=messages,
                                         poisson.llik=TRUE)
  S.sim <- S.sim.res$samples
  n.sim <- dim(S.sim)[1]
  my_shp$pMean_ARR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]-mu0), 1, mean))
  my_shp$pSD_ARR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,]-mu0)), 1, sd)
  my_shp$pMean_RR <- exp(apply(sapply(1:n.sim, function(x) S.sim[x,]), 1, mean))
  my_shp$pSD_RR <- apply(sapply(1:n.sim, function(x) exp(S.sim[x,])), 1, sd)
  out$S.draw <- S.sim
  out$my_shp <- my_shp
  out$para_est <- para_est
  ################################## continuous
  if(is.null(pred.loc)) {
    bound <- raster::aggregate(my_shp)
    regpts <- sp::spsample(bound, cellsize=cellsize, type = "regular")
    vvv <- sp::coordinates(regpts)
    pred.loc <- data.frame(x=vvv[,1], y=vvv[,2])/divisor
    out$bound <- bound
  }
  n.pred.loc <- nrow(pred.loc)
  U.pred <- as.matrix(dist(pred.loc))
  #Sigma x
  Sigma.x2 <- sigma2*exp(-U.pred/phi)
  #########
  #Sigma_x_A
  cov.matrix.x.A=function(pred.loc, S.coord, phi){
    n.pred.loc <- nrow(pred.loc)
    n.distr <- length(S.coord)
    R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
    for (i in 1:n.pred.loc){
      for (j in 1:n.distr){
        U= as.matrix(pdist::pdist(pred.loc[i,],
                                  as.matrix(S.coord[[j]]$xy)))
        R[i,j] =  sum(S.coord[[j]]$weight*exp(-U/phi))
      }
      cat('iter', i, '\n')
    }
    return(R)
  }
  ##########  
  cov.matrix.x.A2=function(pred.loc, S.coord, phi){
    n.pred.loc <- nrow(pred.loc)
    n.distr <- length(S.coord)
    R= matrix(NA, nrow = n.pred.loc, ncol = n.distr)
    for (i in 1:n.pred.loc){
      for (j in 1:n.distr){
        U= as.matrix(pdist::pdist(pred.loc[i,],
                                  as.matrix(S.coord[[j]]$xy)))
        R[i,j] =  mean(exp(-U/phi))
      }
      cat('iter', i, '\n')
    }
    return(R)
  }
  ################
  ##################
  if (weight==TRUE){
    if (parallel==TRUE){
      cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
      NULL
      #Sigma.x.A2 <- sigma2*cov.matrix.x.A4(pred.loc, S.coord, phi)
    }else{
      Sigma.x.A2 <- sigma2*cov.matrix.x.A(pred.loc, S.coord, phi)
    }
  } else{
    if (parallel==TRUE){
      cat("The parallel option will be available once bigstatsr package is submitted on cran, see readme file on github for more info")
      NULL
      #Sigma.x.A2 <- sigma2*cov.matrix.x.A3(pred.loc, S.coord, phi)
    }else{
      Sigma.x.A2 <- sigma2*cov.matrix.x.A2(pred.loc, S.coord, phi)
    }
  }
  #####################
  #Sigma_A
  Sigma.A2 <- Sigma0
  ######
  # matt <- Sigma0
  # diag(matt) <- 1
  # Sigma.A2 <- sigma2*matt
  inv.Sigma.A2 <- solve(Sigma.A2)
  ###############
  ###################
  #The predition
  pred.var2 <- Sigma.x2 - Sigma.x.A2%*%inv.Sigma.A2%*%t(Sigma.x.A2)
  KKK2 <- t(chol(pred.var2))
  S.x2 <- matrix(NA, nrow=n.sim, ncol= n.pred.loc)
  for (i in 1:n.sim){
    pred.mean2 <- Sigma.x.A2%*%(inv.Sigma.A2%*%(S.sim[i,]-mu0))
    S.x2[i,] <- pred.mean2 + KKK2%*%rnorm(n.pred.loc)
  }
  M.E.S.x2 <- exp(apply(S.x2, 2, mean))
  SD.E.S.x2 <- apply(exp(S.x2), 2, sd)
  out$pred.draw <- S.x2
  out$pred.loc <- pred.loc
  out$pred <- M.E.S.x2
  out$predSD <- SD.E.S.x2
  out$call <- match.call()
  attr(out, 'weighted') <- weight
  class(out) <- "SDALGCP"
  return(out)
}
##########################################################################
##' @title Parameter estimation for SDA-LGCP Using Monte Carlo Maximum likelihood
##' @description This function provides the maximum likelihood estimation of the parameter given the precomputed correlation matrices for different values of scale parameter, phi
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param data  data frame containing the variables in the model.
##' @param my_shp A SpatialPolygons orSpatialPolygonsDataFrame  object containing the polygons (i.e each regions).
##' @param delta distance between points
##' @param phi the discretised values of the scale parameter phi. if not supplied, it uses the default.
##' @param pop_shp Optional, The raster of population density map for population weighted approach
##' @param rho Optional, The packing density, default set to 0.55
##' @param weighted To specify if you want to use the population density, default to FALSE, i.e population density is not used.
##' @param plot To display the plot of the points inside the polygon, default to TRUE
##' @param method To specify which method to use to sample the points, the options are 1 for Simple Sequential Inhibition (SSI) process, 2 for Uniform sampling and 3 for regular grid. 1 is the default
##' @param par0 the initial parameter of the fixed effects beta and the variance sigmasq, specified as c(beta, sigma2)
##' @param control.mcmc list from PrevMap package to define the burnin, thining, the number of iteration and the turning parameters see \code{\link{controlmcmcSDA}}.
##' @param plot_profile logical; if TRUE the profile-likelihood is plotted. default is FALSE
##' @details This function performs parameter estimation for a SDA-LGCP Model
##' \bold{Monte Carlo Maximum likelihood.}
##' The Monte Carlo maximum likelihood method uses conditional simulation from the distribution of the random effect \eqn{T(x) = d(x)'\beta+S(x)} given the data \code{y}, in order to approximate the high-dimensiional intractable integral given by the likelihood function. The resulting approximation of the likelihood is then maximized by a numerical optimization algorithm which uses analytic epression for computation of the gradient vector and Hessian matrix. The functions used for numerical optimization are \code{\link{nlminb}} 
##' @return An object of class "PrevMap".
##' The function \code{\link{summary.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{D}: matrix of covariates.
##' @return \code{y}: the count, response observations.
##' @return \code{m}: offset
##' @return \code{beta_opt}: estimates of the fixed effects of the model.
##' @return \code{sigma2_opt}: estimates of the variance of the Gaussian process.
##' @return \code{phi_opt}: estimates of the scale parameter phi of the Gaussian process.
##' @return \code{cov}: covariance matrix of the MCML estimates.
##' @return \code{Sigma_mat_opt}: covariance matrix of the Gaussian process that corresponds to the optimal value
##' @return \code{llike_val_opt}: maximum value of the log-likelihood.
##' @return \code{mu}: mean of the linear predictor
##' @return \code{all_para}: the entire estimates for the different values of phi.
##' @return \code{all_cov}: the entire covariance matrix of the estimates for the different values of phi.
##' @return \code{par0}: the initial parameter of the fixed effects beta and the variance sigmasq used in the estimation
##' @return \code{control.mcmc}: the burnin, thining, the number of iteration and the turning parameters used see \code{\link{controlmcmcSDA}}.
##' @return \code{call}: the matched call.
##' @examples
##' \dontrun{
##' data(PBCshp)
##' data <- as.data.frame(PBCshp@data)
##' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
##' Environment +  offset(log(pop))
##' phi <- seq(500, 1700, length.out = 20)
##' control.mcmc <- list(n.sim = 110000, burnin = 10000, thin= 10, h=1.65/(545^(1/6)),
##'                     c1.h = 0.01, c2.h = 1e-04)
##' my_est <- SDALGCPMCML(formula=FORM, data=data, my_shp=PBCshp, delta=100, phi=phi, method=1, 
##'                      weighted=FALSE,  plot=TRUE, par0=NULL, control.mcmc=control.mcmc)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @importFrom sp bbox
##' @references Giorgi, E., & Diggle, P. J. (2017). PrevMap: an R package for prevalence mapping. Journal of Statistical Software, 78(8), 1-29. doi:10.18637/jss.v078.i08
##' @references Christensen, O. F. (2004). Monte carlo maximum likelihood in model-based geostatistics. Journal of Computational and Graphical Statistics 13, 702-718.
##' @seealso \link{Aggregated_poisson_log_MCML}, \code{\link{Laplace.sampling}}
##' @export
SDALGCPMCML <- function(formula, data, my_shp, delta, phi=NULL, method=1, pop_shp=NULL, 
                        weighted=FALSE,  plot=TRUE, par0=NULL, control.mcmc=NULL, plot_profile=FALSE, rho=NULL){
  if(any(is.na(data))) stop("missing values are not accepted")
  if(class(formula)!="formula") stop("formula must be a 'formula' object that indicates the variables of the fitted model.")
  if(!is.null(control.mcmc) & length(control.mcmc) != 6) stop("please check the input of the controlmcmc argument")
  if(is.null(phi)){
    phi <- seq(sqrt(min(sapply(1:length(my_shp), function(x) my_shp@polygons[[x]]@area))), 
               min(apply(sp::bbox(my_shp), 1, diff))/10, length.out = 20)
  }
  #############create point
  my_list <- SDALGCPpolygonpoints(my_shp=my_shp, delta=delta, method=1, pop_shp=pop_shp, 
                                  weighted=weighted, plot=plot, rho=rho)
  #############precompute matrix
  if(is.null(par0)){
    my_preMatrix <- precomputeCorrMatrix(S.coord = my_list, phi = phi)
  } else{
    phi <- c(phi, par0[length(par0)])
    my_preMatrix <- precomputeCorrMatrix(S.coord = my_list, phi = phi)
  }

  #############estimate parameter
  my_est <- SDAParaEst(formula=formula, data=data, corr= my_preMatrix, par0=par0, 
                       control.mcmc=control.mcmc, plot_profile=plot_profile)
  my_est$call <- match.call()
  attr(my_est, 'SDALGCPMCML') <- TRUE
  class(my_est) <- "SDALGCP"
  return(my_est)
}

##########################################
##' @title Spatial continuous predictions of the relative risk using plug-in of MCML estimates
##' @description This function performs spatial continuous prediction, fixing the model parameters at the Monte Carlo maximum likelihood estimates of a SDA-LGCP model.
##' @param para_est an object of class "SDALGCP" obtained as a result of a call to \code{\link{SDALGCPMCML}}.
##' @param cellsize the size of the computational grid 
##' @param pred.loc optional, the dataframe of the predictive grid.
##' @param continuous logical; to choose which prrediction to do perform, discrete or continuous. the default is continuous.
##' @param control.mcmc output from \code{\link{controlmcmcSDA}}, if not provided, it uses the values used for the parameter estimation
##' @param divisor optional, the value to use to convert the dimension of the polygon, default is 1 which implies no conversion
##' @param plot.correlogram ogical; if plot.correlogram=TRUE the autocorrelation plot of the conditional simulations is displayed.
##' @param messages logical; if messages=TRUE then status messages are printed on the screen (or output device) while the function is running. Default is messages=TRUE.
##' @param parallel to parallelize some part of the function.
##' @details The function returns the prediction of the relative risk exp(S(x))
##' @return pred.draw: the samples of the prediction
##' @return pred: the prediction of the relative risk
##' @return predSD: the standard error of the prediction
##' @return Pred.loc: The cordinates of the predictive locations
##' @examples
##' \dontrun{
##' load("PBCshp")
##' data <- as.data.frame(PBCshp@data)
##' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
##' Environment +  offset(log(pop))
##' phi <- seq(500, 1700, length.out = 20)
##' control.mcmc <- list(n.sim = 110000, burnin = 10000, thin= 10, h=1.65/(545^(1/6)),
##'                     c1.h = 0.01, c2.h = 1e-04)
##' my_est <- SDALGCPMCML(formula=FORM, data=data, my_shp=PBCshp, delta=100, phi=phi, method=1, 
##'                      weighted=FALSE,  plot=TRUE, par0=NULL, control.mcmc=control.mcmc)
##' Con_pred <- SDALGCPPred(para_est=my_est,  cellsize=300, continuous=TRUE)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @references Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). Hierarchical modeling and analysis for spatial data. CRC press.
##' @importFrom sp spsample coordinates
##' @importFrom Matrix solve chol
##' @importFrom pdist pdist 
##' @importFrom stats median
##' @export
SDALGCPPred <- function(para_est, cellsize, continuous=TRUE, control.mcmc=NULL, pred.loc=NULL,
                        divisor=1, plot.correlogram=F, messages=TRUE, parallel=FALSE){
  #############prediction
  if(class(para_est)!="SDALGCP") stop("para_est must be of class 'SDALGCP', that is be an output of SDALGCPMCML function")
  if (continuous==TRUE){
    Con_pred <- SDAContinuousPred(para_est=para_est,  cellsize=cellsize, pred.loc=pred.loc, parallel = parallel, divisor = divisor, 
                                  plot.correlogram = plot.correlogram, messages = messages, control.mcmc = control.mcmc)
  }else{
    Con_pred <- SDADiscretePred(para_est=para_est, control.mcmc = control.mcmc, divisor = divisor, 
                                plot.correlogram = plot.correlogram, messages = messages)
  }
  Con_pred$call <- match.call()
  attr(Con_pred, 'continuous') <- continuous
  class(Con_pred) <- "Pred.SDALGCP"
  return(Con_pred)
}


############################################################################
#######################################
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
###############################
print.SDALGCP <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("Coefficients: \n")
  cf <- c(x$beta_opt, x$sigma2_opt)
  pnames <- pnames <- names(sqrt(diag(x$cov)))
  names(cf) <- pnames
  print(cf)
  cat("\n \n")
  cat("Scale of the spatial correlation, phi: ",x$phi_opt,"\n",sep="")
  cat("Objective function: ",x$llike_val_opt,"\n",sep="")
}

##' @title Summarizing the parameter estimates of SDALGCP model
##' @description \code{summary} method for the class "SDALGCP" that computes the standard errors and p-values of SDALGCP.
##' @param object an object of class "SDALGCP" obatained as result of a call to \code{\link{SDALGCPMCML}} .
##' @param ... further arguments passed to or from other methods.
##' @return A list with the following components
##' @return \code{parameter_estimate_result}: the parameter of the SDALGCP model
##' @return \code{phi}: the scale parameter of the Gaussian process
##' @return \code{ll}: value of likelihood function at the maximum likelihood estimates.
##' @return \code{call}: matched call.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @method summary SDALGCP
##' @export
summary.SDALGCP <- function(object, ...) {
  out <- list()
  cmat  <- cbind(c(object$beta_opt, object$sigma2_opt))
  cmat <- cbind(cmat, sqrt(diag(object$cov)))
  cmat <- cbind(cmat, cmat[, 1]/cmat[, 2])
  parameter_estimate_result <- cbind(cmat, 2*pnorm(-abs(cmat[, 3])))
  colnames(parameter_estimate_result) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  out$parameter_estimate_result  <- parameter_estimate_result
  out$phi <-  object$phi_opt
  out$ll <- object$llike_val_opt
  out$call <- object$call
  class(out) <- "summary.SDALGCP"
  return(out)
}

##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
print.summary.SDALGCP <- function(x, ...){
  cat("Call: \n")
  print(x$call)
  cat("Coefficients: \n")
  printCoefmat(x$parameter_estimate_result, P.values=TRUE,has.Pvalue=TRUE, digits=3)
  cat("\n \n")
  cat("Scale of the spatial correlation, phi: ", x$phi,"\n",sep="")
  cat("Objective function: ", x$ll, "\n",sep="")
  cat("Legend: \n")
  cat("sigma^2 is the variance of the Gaussian process")
}
##########################################
##' @title Plot of the spatially discrete prediction over the region
##' @description This function plots the region specific incidence and covariate adjusted relative risk and their corresponding standard error using the output from \code{\link{SDALGCPPred}} function.
##' @param obj an object of class "Pred.SDALGCP" obtained as result of a call to \code{\link{SDALGCPPred}}
##' @param type Character string: what type of plot to produce. Choices are "incidence" (=exp(mu+S)); "SEincidence" (standard error of incidence); "CovAdjRelRisk" (=exp(S)); or "SECovAdjRelRisk" (standard error of covariate adjusted relative risk);. 
##' @param ... further arguments passed to \code{\link{plot}}.
##' @seealso \code{\link{SDALGCPPred}}
##' @return The function does not return any value.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom sp spplot
##' @export
plot_discrete <- function(obj, type='incidence', ...){
  obj$my_shp$incidence  <- obj$my_shp$pMean_RR
  obj$my_shp$SEincidence  <- obj$my_shp$pSD_RR
  obj$my_shp$CovAdjRelRisk  <- obj$my_shp$pMean_ARR
  obj$my_shp$SECovAdjRelRisk  <- obj$my_shp$pSD_ARR
  sp::spplot(obj$my_shp, type, ...)
}

###################################################################
##' @title Plot of the spatially continuous prediction over the entire region
##' @description This function plots the region specific incidence and covariate adjusted relative risk and their corresponding standard error using the output from \code{\link{SDALGCPPred}} function.
##' @param obj an object of class "Pred.SDALGCP" obtained as result of a call to \code{\link{SDALGCPPred}}
##' @param type Character string: what type of plot to produce. Choices are "relrisk" (=exp(S)); "SErelrisk" (standard error of the relative risk).
##' @param bound the boundary of the predictive grid, not required if predictive grid is not supplied
##' @param ... further arguments passed to \code{\link{plot}}.
##' @seealso \code{\link{SDALGCPPred}}
##' @return The function does not return any value.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom sp spplot
##' @export
plot_continuous <- function(obj, bound=NULL, type='relrisk', ...){
  obj$relrisk  <- obj$pred
  obj$SErelrisk <- obj$predSD
  if (type=="relrisk"){
    dat <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$pred))
    if (any(names(obj) == "bound"))  {
      r1  <- raster::mask(raster::rasterFromXYZ(dat), obj$bound)
      sp::spplot(r1,  ..., sp.layout=obj$bound)
    } else{
      if(is.null(bound)) stop("please supply the boundary of the region")
      r1  <- raster::mask(raster::rasterFromXYZ(dat), bound)
      sp::spplot(r1,  ..., sp.layout=bound)
    }
  }else if (type=='SErelrisk'){
    dat <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$predSD))
    if (any(names(obj) == "bound")){
      r1  <- raster::mask(raster::rasterFromXYZ(dat), obj$bound)
      sp::spplot(r1, ..., sp.layout=obj$bound)
    } else{
      if(is.null(bound)) stop("please supply the boundary of the region")
      r1  <- raster::mask(raster::rasterFromXYZ(dat), bound)
      sp::spplot(r1,  ..., sp.layout=bound)
    }
  }else if (length(type)==2){
    dat1 <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$pred))
    dat2 <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$predSD))
    if (any(names(obj) == "bound")){
      r1  <- raster::mask(raster::rasterFromXYZ(dat1), obj$bound)
      r2  <- raster::mask(raster::rasterFromXYZ(dat2), obj$bound)
      s <- raster::stack(r1, r2)
      #use names.attr = c('Relative Risk', 'Standard Error of Relative Risk') to name the plot
      sp::spplot(s, colorkey=list(space="bottom"), ..., sp.layout=obj$bound)
    } else{
      if(is.null(bound)) stop("please supply the boundary of the region")
      r1  <- raster::mask(raster::rasterFromXYZ(dat1), bound)
      r2  <- raster::mask(raster::rasterFromXYZ(dat2), bound)
      s <- raster::stack(r1, r2)
      sp::spplot(s, colorkey=list(space="bottom"), ..., sp.layout=bound)
    }
  }
}
################################
##' @title Exceedance probability of the relative risk
##' @description Computes the exceedance probability for a given threshold of the spatially continuous relative risk or the region specific relative risk from the object of class "Pred.SDALGCP".
##' @param obj an object of class "Pred.SDALGCP" obtained as result of a call to \code{\link{SDALGCPPred}}.
##' @param thresholds either a vector  of numbers or a vector of single value.
##' @param continuous logical; TRUE for spatially continuous relative risk and FALSE for region specific relative risk. default is TRUE
##' @return A vector or dataframe(for more than one value of the threshold) of the exceedance probability
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
SDALGCPexceedance <- function(obj , thresholds, continuous=TRUE) {
  if (continuous){
    sim <- exp(obj$pred.draw)
    out <- sapply(thresholds, function(x) apply(sim, 2, function(z) mean(z > x)))
  }else{
    sim <- exp(obj$S.draw)
    out <- sapply(thresholds, function(x) apply(sim, 2, function(z) mean(z > x)))
  }
  #colnames(out) <- thresholds
  return(out)
}

##' @title Plot the exceedance probability map
##' @description Plot the exceedance probability map for a given threshold of the spatially continuous relative risk or the region specific relative risk from the object of class "Pred.SDALGCP".
##' @param obj an object of class "Pred.SDALGCP" obtained as result of a call to \code{\link{SDALGCPPred}}.
##' @param thresholds either a vector  of numbers or a vector of single value.
##' @param bound optional; it gives the boundary of the region, only useful when the predictive location is supplied in \link{SDALGCPPred}.
##' @param continuous logical; TRUE for spatially continuous relative risk and FALSE for region specific relative risk. default is TRUE.
##' @param ... further arguments passed to \link{plot}.
##' @return The plot of the exceedance probability map
##' @seealso \link{SDALGCPPred}
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
plot_SDALGCPexceedance <- function(obj, thresholds, bound=NULL, continuous=TRUE, ...){
  if (continuous){
    obj$prob <- SDALGCPexceedance(obj, thresholds=thresholds, continuous=TRUE)
    dat <- data.frame(x=obj$pred.loc[,1], y=obj$pred.loc[,2], z = as.numeric(obj$prob))
    if (any(names(obj) == "bound"))  {
      r1  <- raster::mask(raster::rasterFromXYZ(dat), obj$bound)
      sp::spplot(r1,  ..., sp.layout=obj$bound)
    } else{
      if(is.null(bound)) stop("please supply the boundary of the region")
      r1  <- raster::mask(raster::rasterFromXYZ(dat), bound)
      sp::spplot(r1,  ..., sp.layout=bound)
    }
  }else{
    obj$my_shp$prob <- SDALGCPexceedance(obj, thresholds=thresholds, continuous=FALSE)
    sp::spplot(obj$my_shp, "prob", ...)
  }
}

###############################
##' @title plot.Pred.SDALGCP function
##' @description Simple plotting function for both discrete and conituous prediction from the object of class "Pred.SDALGCP".
##' @param x an object of class "Pred.SDALGCP" obtained as result of a call to \code{\link{SDALGCPPred}}.
##' @param type Character string: what type of plot to produce. For discrete inference choices are "incidence" (=exp(mu+S)); "SEincidence" (standard error of incidence); "CovAdjRelRisk" (=exp(S)); or "SECovAdjRelRisk" (standard error of covariate adjusted relative risk); while for continuous inference, choices are "relrisk" (=exp(S)); "SErelrisk" (standard error of the relative risk).
##' @param continuous logical; TRUE for spatially continuous relative risk and FALSE for region specific relative risk. default is TRUE
##' @param thresholds optional; (only used if you want to plot the exceedance probability) either a vector  of numbers or a vector of single value.
##' @param bound optional; it gives the boundary of the region, only useful when the predictive location is supplied in \link{SDALGCPPred}
##' @param ... further arguments passed to \link{plot}.
##' @seealso \link{SDALGCPPred}, \link{plot_continuous}, \link{plot_discrete}, \link{plot_SDALGCPexceedance}
##' @return The function does not return any value.
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
plot.Pred.SDALGCP <- function(x,  type='relrisk', continuous=TRUE, thresholds=NULL, bound=NULL, ...){
  if(continuous){
    if(is.null(thresholds)){
      plot_continuous(obj=x, bound=bound, type=type, ...)
    } else {
      plot_SDALGCPexceedance(obj=x , thresholds=thresholds, bound=bound, continuous=continuous, ...)
    }
  }else{
    if (is.null(thresholds)){
      if (type=='relrisk') stop("Since you have made a spatially discrete inference, please set type to be one of these four options, choices are 'incidence' (=exp(mu+S)); 'SEincidence' (standard error of incidence); 'CovAdjRelRisk' (=exp(S)); or 'SECovAdjRelRisk' (standard error of covariate adjusted relative risk)")
      plot_discrete(obj=x, type=type, ...)
    } else{
      plot_SDALGCPexceedance(obj=x , thresholds=thresholds, bound=bound, continuous=continuous, ...)
    }
  }
}

##########################
##' @title Confidence Intervals for SDA LGCP Model Parameters
##' @description Computes confidence intervals for one or more parameters in a fitted SDALGCP model from the object of class "SDALGCP", based on asymptotic normality.
##' @param object an object of class "SDALGCP" obtained as result of a call to \code{\link{SDALGCPMCML}}.
##' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
##' @param level the confidence level required.
##' @param dp the number of decimal places for the result
##' @param ... additional argument(s) for methods.
##' @seealso \link{confint.lm}, \link{confint.default}, \link{SDALGCPMCML}
##' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
##' @examples
##' \dontrun{
##' data(PBCshp)
##' data <- as.data.frame(PBCshp@data)
##' FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime + 
##' Environment +  offset(log(pop))
##' phi <- seq(500, 1700, length.out = 20)
##' control.mcmc <- list(n.sim = 110000, burnin = 10000, thin= 10, h=1.65/(545^(1/6)),
##'                     c1.h = 0.01, c2.h = 1e-04)
##' my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=100, phi=phi, method=1, 
##'                      weighted=FALSE,  plot=TRUE, par0=NULL, control.mcmc=control.mcmc)
##' confint(my_est)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
confint.SDALGCP <- function(object, parm, level = 0.95, dp=3, ...){
  cf <- c(object$beta_opt, object$sigma2_opt)
  pnames <- names(sqrt(diag(object$cov)))
  names(cf) <- pnames
  if(missing(parm)) parm <- pnames
  else if(is.numeric(parm)) parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3),
               "%")
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,pct))
  ses <- sqrt(diag(object$cov))[parm]
  ci[] <- cf[parm] + ses %o% fac
  round(ci, dp)
}


##' @title SDA control.mcmc
##' @description This function helps to define the number of iteration, burn-in, thining, and the tunning parameters of the adaptive MALA
##' @param n.sim the number of iteration
##' @param burnin The number of burn-in
##' @param thin the number of thining
##' @param h tuning parameter of the proposal distribution used in the Langevin-Hastings MCMC algorithm (see Laplace.sampling); default is h=NULL and then set internally as 1.65/n(1/6), where n is the dimension of the random effect.
##' @param c1.h value of c1 used in the adaptive scheme for h; default is c1.h=0.01. See also 'Details' in PrevMap package
##' @param c2.h value of c2 used in the adaptive scheme for h; default is c2.h=1e-04. See also 'Details' in PrevMap package
##' @details To be used in SDA
##' @return A list with processed arguments to be passed to the main function.
##' @examples
##' \dontrun{
##' n <- 545
##' h <- 1.65/(n^(1/6))
##' control.mcmc <- controlmcmcSDA(n.sim = 10000, burnin = 2000, 
##' thin= 8, h=h, c1.h = 0.01, c2.h = 1e-04)
##' str(control.mcmc)
##' }
##' @author Olatunji O. Johnson \email{o.johnson@@lancaster.ac.uk}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @importFrom pdist pdist
##' @seealso \link[PrevMap]{control.mcmc.MCML}
##' @export
##' 

controlmcmcSDA <- function(n.sim, burnin, thin, h, c1.h, c2.h){
  control.mcmc <- list(n.sim = n.sim, burnin = burnin, thin= thin, h=h, c1.h = c1.h, c2.h = c2.h)
  return(control.mcmc)
}
