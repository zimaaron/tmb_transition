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
library(colorspace)

## #############################################################################

plot.dir = '/homes/azimmer/comp_plots'

tmb.rd  <- "2018_06_05_14_23_54" 
inla.rd <- "2018_06_05_14_23_53" 

plot_name <- 'HIV - SSSA - GP + NUG'

reg <- 'sssa'

year_list <- 2000:2016

indicator_group <- 'hiv'
indicator <- 'hiv_test'

## now we can compare
tmb.dir   <- sprintf('/share/geospatial/mbg/%s/%s/output/%s', indicator_group, indicator, tmb.rd)
inla.dir  <- sprintf('/share/geospatial/mbg/%s/%s/output/%s', indicator_group, indicator, inla.rd)


## from the simulation, here's the truth param
param.table <- as.data.table(matrix(ncol = 3, nrow = 6))
names(param.table) <- c('param', 'INLA', 'TMB')
param.table$param <- c('int', 'nom. range', 'nom. var', 'time rho', 'nug prec', 'nug sd')

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
inla.params[which(param.table$param == 'nug prec')]   <- summary(inla.model.fit)$hyper[4, 1]
inla.params[which(param.table$param == 'nug sd')]     <- 1 / sqrt(summary(inla.model.fit)$hyper[4, 1])


param.table[, INLA := inla.params]

## grab and convertn TMB estimates
tmb.params <- rep(NA, nrow(param.table))
tmb.params[match(tmb.model.fit$fenames,param.table$param)] <- tmb.model.fit$opt$par[names( tmb.model.fit$opt$par ) == 'alpha_j']
t.lt <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'logtau')] ## tmb log tau
t.lk <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'logkappa')] ## tmb log kappa
t.tr <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'trho')] ## tmb transformed time rho
t.ls <- tmb.model.fit$opt$par[which(names(tmb.model.fit$opt$par) == 'log_nugget_sigma')] ## nug_log_sigma
tmb.params[which(param.table$param == 'nom. range')] <- sqrt(8) / exp(t.lk)
tmb.params[which(param.table$param == 'nom. var')]   <- 1 / (4 * pi * exp(t.lk) ^ 2 * exp(t.lt) ^ 2)
tmb.params[which(param.table$param == 'time rho')] <- (exp(t.tr) - 1) /  (exp(t.tr) + 1)
tmb.params[which(param.table$param == 'nug prec')] <- 1 / exp(t.ls) ^ 2
tmb.params[which(param.table$param == 'nug sd')]  <- exp(t.ls)

param.table[, TMB := tmb.params]

###########################
###########################
#### make some plots ######
###########################
###########################

## load in cell preds
load(paste0(tmb.dir, '/', indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
tmb.cp <- cell_pred
load(paste0(inla.dir, '/', indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
inla.cp <- cell_pred
rm(cell_pred);for(i in 1:10) gc()

## Load simple polygon template to model over
gaul_list           <- get_gaul_codes(reg)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

## define summary stats to calc across cell draws and make summary rasters
summstats <- c('mean', 'sd', 'lower', 'upper', 'cirange')
ss.list <- as.list(summstats)
names(ss.list) <- summstats

for(type in c('inla', 'tmb')){
  message(sprintf('SUMMARIZING %s', type))
  assign(paste0(type, '_raster_list'),
         lapply(ss.list,
                function(x){
                  message(sprintf('-- now finding %s', x))
                  make_cell_pred_summary(get(paste0(type, '.cp')),
                                         mask = simple_raster,
                                         summary_stat = x)
                }
                )
         )
}

## also load in input data (should be same from INLA and TMB)
df <- fread(paste0(inla.dir, '/input_data_bin0_', reg, '_0.csv'))


###################################
## plot mean and cirange rasters ##
###################################

tmb.mean.rast <- tmb_raster_list[['mean']]
tmb.sd.rast   <- tmb_raster_list[['sd']]

inla.mean.rast <- inla_raster_list[['mean']]
inla.sd.rast   <- inla_raster_list[['sd']]

## png(sprintf('%s/inla_tmb%s_%s_%s_plot.png',plot.dir, plot_name, inla.rd, tmb.rd), height=60,width=18, units = 'in', res = 350)

pdf(sprintf('%s/inla_tmb%s_%s_%s_plot.pdf',plot.dir, plot_name, inla.rd, tmb.rd), height=60,width=18)

grid.table(param.table)

for(thing in c('mean','sd')){

  
    layout(matrix(1:(17 * 5), 17, 5, byrow = TRUE))
    for(i in 1:17){

      message(sprintf('On year %i', year_list[i]))

      ## select some random pixels to scatter
      samp=sample(cellIdx(inla.mean.rast[[1]]),1e4)

      ## subset data to 1 year
      tmp <- subset(df, year == year_list[i])
      
      ## grab this years summaries from inla and tmb models
      rinla <-  inla_raster_list[[thing]][[i]]
      rtmb  <-  tmb_raster_list[[thing]][[i]]

      ## plot each raster, their difference, and a scatter of a subset of pixels
      par(mar = c(0, 0, 1.4, 1),bty='n')
      maxes <- max(c(as.vector(rtmb),as.vector(rinla)),na.rm=TRUE)

      ## summary rasters
      plot(rtmb,  maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend.width = 5, legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main='TMB',zlim=c(0,maxes))
      plot(rinla, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=FALSE, main='R-INLA',zlim=c(0,maxes))
      
      ## difference
      plot(rinla-rtmb, maxpixel=1e7, col=rainbow_hcl(100), axes=FALSE, legend=T, legend.width = 5, main=paste0('DIFFERENCE (I-T): ',thing))
      
      ## scatter
      par(mar = c(4, 4, 2, 2),bty='n')
      plot(x=as.vector(rinla)[samp],y=as.vector(rtmb)[samp],
           xlab='R-INLA',ylab='TMB',cex=.01,pch=19,main='Pixel Scatter',xlim = c(0, maxes), ylim = c(0, maxes))
      lines(x=c(0,maxes),y=c(0,maxes),col='red')

      ## plot data points over the shapefile
      plot(simple_polygon, main='Data Locs')
      if(nrow(tmp) > 0){
        points(x=tmp$long,y=tmp$lat, pch=19, cex=(tmp$N / max(tmp$N)) )
      }
      
    }
}

dev.off()


######################
## plot time series ##
######################

t.agg0 <- fread(sprintf('%s/pred_derivatives/admin_summaries/%s_admin_0_unraked_summary.csv', tmb.dir, indicator))
t.agg0[, level := 0]
t.agg1 <- fread(sprintf('%s/pred_derivatives/admin_summaries/%s_admin_1_unraked_summary.csv', tmb.dir, indicator))
t.agg1[, level := 1]
t.agg2 <- fread(sprintf('%s/pred_derivatives/admin_summaries/%s_admin_2_unraked_summary.csv', tmb.dir, indicator))
t.agg2[, level := 2]
t.agg  <- as.data.table(rbind.fill(t.agg0, t.agg1, t.agg2))
t.agg[, model := 'TMB']

i.agg0 <- fread(sprintf('%s/pred_derivatives/admin_summaries/%s_admin_0_unraked_summary.csv', inla.dir, indicator))
i.agg0[, level := 0]
i.agg1 <- fread(sprintf('%s/pred_derivatives/admin_summaries/%s_admin_1_unraked_summary.csv', inla.dir, indicator))
i.agg1[, level := 1]
i.agg2 <- fread(sprintf('%s/pred_derivatives/admin_summaries/%s_admin_2_unraked_summary.csv', inla.dir, indicator))
i.agg2[, level := 2]
i.list <- list(a0 = i.agg0,
               a1 = i.agg1,
               a2 = i.agg2)
i.agg  <- as.data.table(rbind.fill(i.agg0, i.agg1, i.agg2))
i.agg[, model := 'INLA']

all.agg <- rbind(t.agg, i.agg)

## now we also need to get the true time series to add to the plot


pdf(sprintf('%s/inla_tmb_%s_time_aggregation_comparison_%s_%s_plot.pdf',plot.dir, plot_name, inla.rd, tmb.rd), height=12,width=18)

for(ll in 0:1){

  plot.dat <- subset(all.agg, level == ll)
  p <- ggplot(plot.dat) + geom_line( aes(x = year, y = mean, group = model, color = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = year, group = model, fill = model), alpha = 0.2)
  

  if(nrow(unique(plot.dat[, paste0('ADM', ll, '_CODE'), with = F])) > 1){
    p <- p + facet_wrap( as.formula(paste0('~ ADM', ll, '_NAME') ))
  }else{
    p <- p + ggtitle(unique(plot.dat[, paste0('ADM', ll, '_NAME'), with = F]))
  }
  plot(p)
  
}

dev.off()

