## some functions to help with realistic inla/tmb simulation comparison
## written by a0z 5/17/18


## first, a function to simulate from a GF using INLA
## this function is taken from the spde tutorial
rspde <- function (coords, kappa, variance = 1, alpha = 2, n = 1, mesh,
                   verbose = FALSE, seed, return.attributes = FALSE)
{
    t0 <- Sys.time()
    theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
    if (verbose)
        cat("theta =", theta, "\n")
    if (missing(mesh)) {
        mesh.pars <- c(0.5, 1, 0.1, 0.5, 1) * sqrt(alpha - ncol(coords)/2)/kappa
        if (verbose)
            cat("mesh.pars =", mesh.pars, "\n")
        attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords),
                                                        ],
                                               max.edge = mesh.pars[1:2],
                                               cutoff = mesh.pars[3],
                                               offset = mesh.pars[4:5]))
        if (verbose)
            cat("n.mesh =", attributes$mesh$n, "\n")
    }
    else attributes <- list(mesh = mesh)
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
    attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
    attributes$A <- inla.mesh.project(mesh = attributes$mesh,
                                      loc = coords)$A
    if (n == 1)
        result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q,
                                                     constr = attributes$spde$f$extraconstr))
    t1 <- Sys.time()
    result <- inla.qsample(n, attributes$Q, seed = ifelse(missing(seed),
                                                          0, seed),
                           constr = attributes$spde$f$extraconstr)
    if (nrow(result) < nrow(attributes$A)) {
        result <- rbind(result, matrix(NA, nrow(attributes$A) -
                                           nrow(result), ncol(result)))
        dimnames(result)[[1]] <- paste("x", 1:nrow(result), sep = "")
        for (j in 1:ncol(result)) result[, j] <- drop(attributes$A %*%
                                                      result[1:ncol(attributes$A),
                                                             j])
    }
    else {
        for (j in 1:ncol(result)) result[1:nrow(attributes$A),
                                         j] <- drop(attributes$A %*% result[, j])
        result <- result[1:nrow(attributes$A), ]
    }
    t2 <- Sys.time()
    attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 -
                                                              t0)
    if (return.attributes)
        attributes(result) <- c(attributes(result), attributes)
    return(drop(result))
}


sim.realistic.data <- function(reg,
                               year_list,
                               betas,
                               sp.kappa,
                               sp.alpha,
                               t.rho,
                               n.clust,
                               m.clust,
                               covs,
                               simple_raster,
                               simple_polygon, 
                               out.dir){

########################################
## load and prepare covariate rasters ##
########################################  
  
  cov_layers <- load_and_crop_covariates_annual(covs            = covs[, name],
                                                measures        = covs[, meas],
                                                simple_polygon  = simple_polygon,
                                                start_year      = min(year_list),
                                                end_year        = max(year_list),
                                                interval_mo     = 12) ## always grab annual, then subset if need be
  ## subset to years in yearlist
  for(cc in 1:length(cov_layers)){
    if(dim(cov_layers[[cc]])[3] > 1){
      cov_layers[[cc]] <- cov_layers[[cc]][[which( min(year_list):max(year_list) %in% year_list )]]
    }

    ## center-scale covs. (TODO by year or across years?)
    cov_layers[[cc]] <- (cov_layers[[cc]] - mean(values(cov_layers[[cc]]), na.rm = T)) / sd(values(cov_layers[[cc]]), na.rm = T)
  }

  ## we also want our cov_layers to align with simple_raster
  for(l in 1:length(cov_layers)) {
    message(sprintf("On cov %i out of %i", l, length(cov_layers)))
    cov_layers[[l]]  <- crop(cov_layers[[l]], extent(simple_raster))
    cov_layers[[l]]  <- setExtent(cov_layers[[l]], simple_raster)
    cov_layers[[l]]  <- mask(cov_layers[[l]], simple_raster)
  }


  ## plot center-scaled covariates
  dir.create(sprintf('%s/simulated_obj/', out.dir))
  pdf(sprintf('%s/simulated_obj/cov_plot.pdf', out.dir), width = 16, height = 16)
  for(cc in 1:length(cov_layers)){
    message(sprintf('Plotting covariate: %s\n', names(cov_layers)[[cc]]))
    if(dim(cov_layers[[cc]])[3] == 1){
      message('--plotting single synoptic map\n')
      par(mfrow = c(1, 1))
    }else{
      message('--plotting time series\n')
      par(mfrow = rep( ceiling(sqrt( dim(cov_layers[[cc]])[3] )), 2))
    }
    for(yy in 1:dim(cov_layers[[cc]])[3]){
      raster::plot(cov_layers[[cc]][[yy]],
                   main = paste(names(cov_layers)[cc],
                                ifelse(dim(cov_layers[[cc]])[3] == 1,
                                       'synoptic',
                                       year_list[yy]), 
                                sep = ': '))
    }
  }
  dev.off()

  ## now we can simulate our true surface

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ############################
  ## simulate space-time gp ##
  ############################

  ## to simulate, we need lat-lon locs for the entire raster

  ## convert simple raster of our region to spatialpolygonsDF
  pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)

  ## reproject sp obj
  geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
  pix.pts <- spTransform(pix.pts, CRS(geo.prj)) 
  proj4string(pix.pts)

  ## get coords
  pix.pts@data <- data.frame(pix.pts@data, long=coordinates(pix.pts)[,1],
                             lat=coordinates(pix.pts)[,2])
  pix.pts.numeric <- as.data.frame(pix.pts@data)

  ## now we can use these coords to simulate GP from rspde()
  reg.mesh <- inla.mesh.2d(boundary = inla.sp2segment(simple_polygon),
                           loc = pix.pts@data[, 2:3],
                           max.edge = c(0.25, 5),
                           offset = c(1, 5),
                           cutoff = 0.25)


  ## get gp fields across space that are indep in time
  gp.iid <- rspde(coords = cbind(pix.pts.numeric[, 2], pix.pts.numeric[, 3]),
                  kappa = sp.kappa,
                  variance = sp.var,
                  alpha = sp.alpha,
                  mesh = reg.mesh,
                  n = length(year_list))

  ## introduce temporal ar1 correlation at the pixel level
  gp.cor <- gp.iid
  for(ii in 2:ncol(gp.cor)){
    gp.cor[, ii] <- t.rho * gp.cor[, ii - 1] + sqrt(1 - t.rho ^ 2) * gp.iid[, ii]
  }

  
  ## convert them to rasters
  for(cc in 1:ncol(gp.iid)){
    if(cc == 1){
      gp.rast <- rasterize(x = pix.pts@data[, 2:3],
                           y = simple_raster,
                           field = gp.cor[, cc])
    }else{
      gp.rast <- addLayer(gp.rast,
                          rasterize(x = pix.pts@data[, 2:3],
                                    y = simple_raster,
                                    field = gp.cor[, cc]))
    }
  }

  pdf(sprintf('%s/simulated_obj/st_gp_plot.pdf', out.dir), width = 16, height = 16)
  par(mfrow = rep( ceiling(sqrt( dim(gp.rast)[3] )), 2))
  for(yy in 1:dim(gp.rast)[3]){
    raster::plot(gp.rast[[yy]],
                 main = paste('GP',
                              year_list[yy]), 
                 sep = ': ')
  }
  dev.off()

  ## now we can make the true surface

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #######################################################
  ## make true surface by combining cov effects and gp ##
  #######################################################

  ## finally, we combine the gp and the covariate effecgts to get our surface in link (e.g. logit if binomial) space
  true.rast <- gp.rast
  for(cc in 1:length(cov_layers)){
    true.rast <- true.rast + betas[cc] * cov_layers[[cc]] ## should work for both stationary and time-varying
  }

  ## we append the gp to the cov_layers
  cov_layers[['gp']] <- gp.rast


  pdf(sprintf('%s/simulated_obj/true_surface_plot.pdf', out.dir), width = 16, height = 16)
  par(mfrow = rep( ceiling(sqrt( dim(true.rast)[3] )), 2))
  for(yy in 1:dim(true.rast)[3]){
    raster::plot(true.rast[[yy]],
                 main = paste('GP',
                              year_list[yy]), 
                 sep = ': ')
  }
  dev.off()

  ## now the surface simulation is done and all we need to do is simulate the data

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #####################################
  ## simulate data from true surface ##
  #####################################

  ## randomly (for now) select data boservation locations across time

  ## to do this, we sample, with replacement, from the lat-longs that we used to sim the GP
  sim.rows <- sample(x = 1:nrow(pix.pts.numeric), size = n.clust * length(year_list),
                     replace = TRUE)
  sim.dat <- as.data.table(pix.pts.numeric[, -1])
  sim.dat <- sim.dat[sim.rows, ]

  ## add in years
  sim.dat[, year := rep(year_list, each = n.clust)]

  ## extract the value of the true surface at data locations
  true_p_logit<- numeric(nrow(sim.dat))
  for(yy in unique(year_list)){
    true_p_logit[which(sim.dat[, year] == yy)] <- raster::extract(x = true.rast, y = sim.dat[year == yy, .(long, lat)],
                                                                  layer = which(year_list %in% yy))
  }
  sim.dat[, p_true := inv.logit(true_p_logit)]

  ## add in cluster sample size
  sim.dat[, N := rpois(n = nrow(sim.dat), lambda = m.clust)]

  ## and now we simulate binomial observations from the true surface
  sim.dat[, Y := rbinom(n = nrow(sim.dat), size = sim.dat[, N], prob = sim.dat[, p_true])]

  ## and get empirical p_obs
  sim.dat[, p_obs := Y / N]

  ## now we just finish by making some convenience objects and saving everything

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ############################################################################################
  ## for convenience, we also extract covariate values to the same df (and the true gp val) ##
  ############################################################################################
  cov.mat <- matrix(ncol = length(cov_layers),
                    nrow = nrow(sim.dat))
  for( cc in 1:length(cov_layers) ){
    tmp <- numeric(nrow(sim.dat))

    if(dim(cov_layers[[cc]])[3] > 1){
      for(ll in 1:dim(cov_layers[[cc]])[3]){
        tmp[which(sim.dat[, year] == year_list[ll])] <- raster::extract(x = cov_layers[[cc]], y = sim.dat[year == year_list[ll], .(long, lat)],
                                                                        layer = ll)[, 1]
      }
    } else{
      tmp <- raster::extract(x = cov_layers[[cc]], y = sim.dat[, .(long, lat)])
    }

    cov.mat[, cc] <- tmp  
  }
  cov.mat <- as.data.table(cov.mat)
  setnames(cov.mat, names(cov_layers))

  #####################################
  ## combine into a single master df ##
  #####################################
  sim.dat <- cbind(sim.dat, cov.mat)


  ###################################
  ## save everything we might want ##
  ###################################
  saveRDS(object = sim.dat,
          file = sprintf('%s/simulated_obj/sim_data.rds', out.dir))

  saveRDS(object = cov_layers,
          file = sprintf('%s/simulated_obj/cov_gp_rasters.rds', out.dir))

  saveRDS(object = reg.mesh,
          file = sprintf('%s/simulated_obj/region_mesh.rds', out.dir))

  #########################
  ## return a named list ##
  #########################
  return(list(sim.dat = sim.dat,
              cov_gp_raster = cov_layers, 
              mesh_s = reg.mesh))
}
