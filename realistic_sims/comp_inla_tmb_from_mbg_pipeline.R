## this script compares two simulated model outputs from our full mbg pipeline.
## it is in this dir since I'm using it to compare differences between model fits in INLA and TMB
## aoz

## #############################################################################
## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_transition', user)
pull_tmb_git <- FALSE

## grab libraries and functions from MBG code
setwd(core_repo)
commondir    <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, 'mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)


## load packages
require(raster)
require(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

## #############################################################################

plot.dir = '/homes/azimmer'

tmb.rd  <- "2018_05_22_10_45_31"
## inla.rd <- "2018_05_22_10_45_29" ## no nug, data all years, all year time mesh
inla.rd <- '2018_05_23_15_34_59' ## no nug, data all years, c(1, 6, 11, 16) time mesh
true.rd <- "2018_05_18_13_03_51"

reg = 'nga'

indicator <- indicator_group <- 'simulation'


## now we can compare
tmb.dir   <- sprintf('/share/geospatial/mbg/%s/%s/output/%s', indicator_group, indicator, tmb.rd)
inla.dir  <- sprintf('/share/geospatial/mbg/%s/%s/output/%s', indicator_group, indicator, inla.rd)
truth.dir <- sprintf('/homes/azimmer/tmb_inla_sim/%s/simulated_obj', true.rd)


## from the simulation, here's the truth params
param.table <- readRDS(paste0(truth.dir, '/true_param_table.rds'))

## load in some model run stuff from INLA and TMB
inla.model.fit <- readRDS( file = sprintf('/share/geospatial/mbg/%s/%s/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                          indicator_group, indicator, inla.rd, 'inla', reg, 0, 0))
tmb.model.fit  <- readRDS( file = sprintf('/share/geospatial/mbg/%s/%s/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                          indicator_group, indicator, tmb.rd, 'tmb', reg, 0, 0))

## grab and convert INLA estimates
inla.params <- rep(NA, nrow(param.table))
inla.params[match(rownames(summary(inla.model.fit)$fixed), param.table$param)] <- summary(inla.model.fit)$fixed[, 1]
i.lt <- summary(inla.model.fit)$hyper[1, 1] ## logtau
i.lk <- summary(inla.model.fit)$hyper[2, 1] ## logkappa
inla.params[which(param.table$param == 'nom. range')] <- sqrt(8) / exp(i.lk)
inla.params[which(param.table$param == 'nom. var')]   <- 1 / (4 * pi * exp(i.lk) ^ 2 * exp(i.lt) ^ 2)
inla.params[which(param.table$param == 'time rho')]   <- summary(inla.model.fit)$hyper[3, 1]

param.table[, INLA := inla.params]

## grab and convertn TMB estimates
tmb.params <- rep(NA, nrow(param.table))
tmb.params[match(tmb.model.fit$fenames,param.table$param)] <- tmb.model.fit$opt$par[names( tmb.model.fit$opt$par ) == 'alpha_j']
t.lt <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'logtau')] ## tmb log tau
t.lk <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'logkappa')] ## tmb log kappa
t.tr <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'trho')] ## tmb transformed time rho
tmb.params[which(param.table$param == 'nom. range')] <- sqrt(8) / exp(t.lk)
tmb.params[which(param.table$param == 'nom. var')]   <- 1 / (4 * pi * exp(t.lk) ^ 2 * exp(t.lt) ^ 2)
tmb.params[which(param.table$param == 'time rho')] <- (exp(t.tr) - 1) /  (exp(t.tr) + 1)

param.table[, TMB := tmb.params]

###########################
###########################
#### make some plots ######
###########################
###########################

###################################
## plot mean and cirange rasters ##
###################################

tmb.mean.rast <- brick(paste0(tmb.dir, '/simulation_nga_mean_raster.tif'))
tmb.cirange.rast <- brick(paste0(tmb.dir, '/simulation_nga_cirange_raster.tif'))
inla.mean.rast <- brick(paste0(inla.dir, '/simulation_nga_mean_raster.tif'))
inla.cirange.rast <- brick(paste0(inla.dir, '/simulation_nga_cirange_raster.tif'))
true.rast <- readRDS(paste0(truth.dir, '/true_surface_raster.rds'))
dt <- readRDS(paste0(truth.dir, '/sim_data.rds'))
dt$period_id = as.numeric(as.factor(dt$year))


## load in the region shapefile and prep the boundary
gaul_list           <- get_gaul_codes(reg)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]


pdf(sprintf('%s/inla_tmb_full_pipeline_raster_comparison_%s_%s_plot.pdf',plot.dir, inla.rd, tmb.rd), height=15,width=30)

grid.table(param.table)

plot.in.logit.space <- FALSE

for(thing in c('mean','cirange')){
    layout(matrix(1:24, 4, 6, byrow = TRUE))
    samp=sample(cellIdx(inla.mean.rast[[1]]),1e4)
    for(i in seq(1, 16, by = 5)){

      ## TODO: allow plotting in logit space
      if(plot.in.logit.space){
        
      }else{
        
        if(thing=='mean'){
          rinla <-  inla.mean.rast[[i]]
          rtmb  <-  tmb.mean.rast[[i]]
        }
        if(thing=='cirange'){
          rinla <-  inla.cirange.rast[[i]]
          rtmb  <-  tmb.cirange.rast[[i]]
        }

        true <- true.rast[[i]] ## in logit space
        values(true) <- plogis(values(true))
      }

      
      tmp <- subset(dt, period_id==i)
      
      # rasters
      par(mar = c(0, 0, 1.4, 1),bty='n')
      maxes <- max(c(as.vector(rtmb),as.vector(rinla)),na.rm=TRUE)
      plot(rinla-rtmb, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=T, main=paste0('DIFFERENCE ',thing))
      plot(rtmb,  maxpixel=1e7, col=rainbow(100), axes=FALSE, legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main='TMB',zlim=c(0,maxes))
      plot(rinla, maxpixel=1e7, col=rainbow(100), axes=FALSE, legend=FALSE, main='R-INLA',zlim=c(0,maxes))
      plot(true, maxpixel=1e7, col=rainbow(100), axes=FALSE, legend=FALSE, main='TRUE (mean)',zlim=c(0,maxes))
      
      # scatter
      par(mar = c(4, 4, 2, 2),bty='n')
      plot(x=as.vector(rinla)[samp],y=as.vector(rtmb)[samp],xlab='R-INLA',ylab='TMB',cex=.01,pch=19,main='COMPARE')
      lines(x=c(0,maxes),y=c(0,maxes),col='red')

      # plot data points over the shapefile
      plot(simple_polygon, main='DATA LOCATIONS')
      points( x=tmp$long,y=tmp$lat, pch=19, cex=(tmp$N / max(tmp$N)) )
      
    }
}

dev.off()


######################
## plot time series ##
######################

t.agg0 <- fread(sprintf('%s/pred_derivatives/admin_summaries/simulation_admin_0_unraked_summary.csv', tmb.dir))
t.agg0[, level := 0]
t.agg1 <- fread(sprintf('%s/pred_derivatives/admin_summaries/simulation_admin_1_unraked_summary.csv', tmb.dir))
t.agg1[, level := 1]
t.agg2 <- fread(sprintf('%s/pred_derivatives/admin_summaries/simulation_admin_2_unraked_summary.csv', tmb.dir))
t.agg2[, level := 2]
t.agg  <- as.data.table(rbind.fill(t.agg0, t.agg1, t.agg2))
t.agg[, model := 'TMB']

i.agg0 <- fread(sprintf('%s/pred_derivatives/admin_summaries/simulation_admin_0_unraked_summary.csv', inla.dir))
i.agg0[, level := 0]
i.agg1 <- fread(sprintf('%s/pred_derivatives/admin_summaries/simulation_admin_1_unraked_summary.csv', inla.dir))
i.agg1[, level := 1]
i.agg2 <- fread(sprintf('%s/pred_derivatives/admin_summaries/simulation_admin_2_unraked_summary.csv', inla.dir))
i.agg2[, level := 2]
i.list <- list(a0 = i.agg0,
               a1 = i.agg1,
               a2 = i.agg2)
i.agg  <- as.data.table(rbind.fill(i.agg0, i.agg1, i.agg2))
i.agg[, model := 'INLA']

all.agg <- rbind(t.agg, i.agg)


pdf(sprintf('%s/inla_tmb_full_pipeline_time_aggregation_comparison_%s_%s_plot.pdf',plot.dir, inla.rd, tmb.rd), height=12,width=18)

for(ll in 0:2){

  plot.dat <- subset(all.agg, level == ll)
  p <- ggplot(plot.dat) + geom_line( aes(x = year, y = mean, group = model, color = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = year, group = model, fill = model), alpha = 0.2)
  

  if(nrow(unique(plot.dat[, paste0('ADM', ll, '_CODE'), with = F])) > 1){
    p <- p + facet_wrap( as.formula(paste0('~ ADM', ll, '_CODE') ))
  }else{
    p <- p + ggtitle(unique(plot.dat[, paste0('ADM', ll, '_CODE'), with = F]))
  }
  plot(p)
  
}

dev.off()

