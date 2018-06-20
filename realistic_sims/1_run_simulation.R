## this script simulates some realistic datasets for comparison between INLA and TMB
## it leverages existing architecture that the LBD team at IHME has already created
## written by AOZ
## last editted 5/14/18

## options(error = recover)
load('/homes/azimmer/scratch/tmb_space_debug.RData')

## DO THIS!
################################################################################
## ADD A NOTE! to help identify what you were doing with this run
logging_note <- 'Running these models with real HIV SSSA data. using weights'
################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par.iter <- 1  ## as.numeric(  commandArgs()[4]) ## all we need is to grab the (parallel) iteration of this run
run_date <- "2018_06_12_15_54_38" ## as.character(commandArgs()[5]) ## and the run_date so we know where to load from

#############################################
## setup the environment for singularity R ##
#############################################

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

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

## Now we can switch to the TMB repo
setwd(tmb_repo)
if(pull_tmb_git) system(sprintf('cd %s\ngit pull %s %s', core_repo, remote, branch))
source('./realistic_sims/realistic_sim_utils.R')

## setup is now done. setup some parameters for this simulation

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################################################
## load in the loopvars from launch and setup all the params for this job ##
############################################################################
out.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', run_date)
loopvars <- readRDS(file = paste0(out.dir, '/loopvars.rds'))

reg             <- as.character(loopvars[par.iter, 1])
year_list       <- eval(parse(text = loopvars[par.iter, 2]))
cov_names       <- eval(parse(text = as.character(loopvars[par.iter, 3])))
cov_measures    <- eval(parse(text = as.character(loopvars[par.iter, 4])))
betas           <- eval(parse(text = as.character(loopvars[par.iter, 5])))
alpha           <- as.numeric(loopvars[par.iter, 6])
sp.range        <- as.numeric(loopvars[par.iter, 7])
sp.var          <- as.numeric(loopvars[par.iter, 8])
sp.alpha        <- as.numeric(loopvars[par.iter, 9])
nug.var         <- as.numeric(loopvars[par.iter, 10])
t.rho           <- as.numeric(loopvars[par.iter, 11])
mesh_s_max_edge <- as.character(loopvars[par.iter, 12])
n.clust         <- as.numeric(loopvars[par.iter, 13])
m.clust         <- as.numeric(loopvars[par.iter, 14])
sample.strat    <- as.character(loopvars[par.iter, 15])
cores           <- as.numeric(loopvars[par.iter, 16]) 
ndraws          <- as.numeric(loopvars[par.iter, 17])

covs <- data.table(name = cov_names, meas = cov_measures)

## TODO? add in some validation options?

## and I hardcode a few other options that are useful sometimes when running interactively

## to make it easier to run real data from this code, usually have SIM=TRUE, REAL=FALSE
use_real_data <- FALSE
use_sim_data  <- TRUE

## should we save this data to mbg input_data? useful if we want to run sim data through lbd_core
save.as.input.data <- FALSE
data.tag <- '_allyrs_nug'

## for use in INLA
ncores <- 3

## end of user inputs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## transform some inputs into other useful quantities

## convert sp.range and sp.var into sp.kappa for rspde function
sp.kappa <- sqrt(8) / sp.range


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################
## load in region/counry shapes and covs ##
###########################################
dir.create(sprintf('%s/simulated_obj', out.dir), recursive = TRUE)

## load in the region shapefile and prep the boundary
gaul_list           <- get_gaul_codes(reg)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

###################
## simulate data ##
###################

## make an object with true param values
true.params <- data.table(param = c('int',
                                    cov_names,
                                    'nom. range',
                                    'nom. var',
                                    'time rho'
                                    ),
                          truth = c(alpha,
                                    betas,
                                    sp.range,
                                    sp.var,
                                    t.rho))

saveRDS(file = sprintf('%s/simulated_obj/true_param_table.rds', out.dir),
        object = true.params)

sim.obj <- sim.realistic.data(reg = reg,
                              year_list = year_list,
                              betas = betas,
                              sp.kappa = sp.kappa,
                              sp.alpha = sp.alpha,
                              t.rho = t.rho,
                              nug.var = nug.var, 
                              n.clust = n.clust,
                              m.clust = m.clust,
                              covs = covs,
                              simple_raster = simple_raster,
                              simple_polygon = simple_polygon, 
                              out.dir = out.dir,
                              seed = 123456)

saveRDS(file = sprintf('%s/simulated_obj/sim_obj.rds', out.dir),
        object = sim.obj)

dt <- sim.obj$sim.dat ## simulated data, lat-long, year, covs, true surface
covs.gp <- sim.obj$cov.gp.rasters   ## rasters of covs and true simulated gp field
true.gp <- covs.gp[['gp']]
cov_list <- covs.gp$gp <- NULL
true.rast <- sim.obj$true.rast

## save (if desired) this simulated dataset to .../mbg/input_data for mbg pipeline
if(save.as.input.data){
  df <- data.table(longitude = dt$long,
                   latitude = dt$lat,
                   year = dt$year,
                   country = reg,
                   N = dt$N,
                   simulation = dt$Y, 
                   weight = 1)
  write.csv(file = sprintf('/share/geospatial/mbg/input_data/simulation%s.csv', data.tag),
            x = df)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################
#####################################
## setup for tmb and INLA modeling ##
#####################################
#####################################
dir.create(sprintf('%s/modeling/inputs', out.dir), recursive = TRUE)
dir.create(sprintf('%s/modeling/tmb/outputs', out.dir), recursive = TRUE)
dir.create(sprintf('%s/modeling/inla/outputs', out.dir), recursive = TRUE)

#########################################
## load in some real data (if desired) ##
#########################################

if(use_real_data){

  reg <- 'sssa'
  indicator = 'hiv_test'
  indicator_group = 'hiv'
  age <- holdout <- test <- 0
  yearload <- 'annual'
  withtag <- TRUE
  datatag <- '_survey'
  use_share <- 'FALSE'
  year_list <- 2000:2016
  
  pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

  
  ## load in the region shapefile and prep the boundary
  gaul_list           <- get_gaul_codes(reg)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]

  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape)
  simple_raster      <- raster_list[['simple_raster']]
  pop_raster         <- raster_list[['pop_raster']]

  run_date <- make_time_stamp(TRUE)
  dt <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                        simple      = simple_polygon,
                        agebin      = age,
                        removeyemen = TRUE,
                        pathaddin   = pathaddin,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        yl          = year_list)

  ## just to make sure everything goes smoothly, add in single datapoints to missing years
  missing.yrs <- setdiff(year_list, unique(dt[, year]))

  if(length(missing.yrs) > 0){
    for(yy in missing.yrs){
      new.row <- dt[1, ]
      new.row[, weight := 0]
      new.row[, year := yy]
      dt <- rbind(dt, new.row)
    }
  }
  

  dt[, Y:= hiv_test]

}


###############################
## SETUP SOME SHARED OBJECTS ##
###############################

## ~~~
## required for model fit
## ~~~

dt[, id := 1:.N]
dt[, period_id := as.numeric(as.factor(dt[, year]))]
dt.coords <- cbind(dt$long, dt$lat)
dt.pers   <- dt[, period_id]
nperiods  <- length(year_list)

## Build spatial mesh over modeling area
mesh_s <- build_space_mesh(d           = dt,
                           simple      = simple_polygon,
                           max_edge    = mesh_s_max_edge,
                           mesh_offset = "c(1, 5)")

pdf(sprintf('%s/modeling/inputs/mesh.pdf', out.dir))
plot(mesh_s)
plot(simple_raster, add = TRUE) ## just to show loc of simple_raster under mesh for scale
plot(mesh_s, add = TRUE)
dev.off()

nodes <- mesh_s$n ## get number of mesh nodes
spde <- inla.spde2.matern( mesh_s,alpha=2 )
## Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
## ^ this gives us a linear reduction of \Sigma^{-1} as:
## \Sigma^{-1} = \kappa^4 M_0 + 2\kappa^2M_1 + M_2
## M_2 = M_1M_0^{-1}M_1
## Where the Ms are all sparse matrices stored as "dgTMatrix"
## names(spde$param.inla)

## use inla helper functions to project the spatial effect from mesh points to data points
A.proj <- inla.spde.make.A(mesh  = mesh_s,
                           loc   = dt.coords,
                           group = dt.pers)

## save relevant objects
saveRDS(file = sprintf('%s/modeling/inputs/mesh.rds', out.dir), mesh_s)
saveRDS(file = sprintf('%s/modeling/inputs/spde.rds', out.dir), spde)

## ~~~
## required for predict
## ~~~

## pull out covariates in format we expect them
## a list of length periods with a brick of named covariates inside
cov_list <- covs.gp
cov_list$gp <- NULL
new_cl <- list()
for(p in 1:nperiods){
  new_cl[[p]] <- list()
  for(n in names(cov_list)){
    if(dim(cov_list[[n]])[3] == 1){
      new_cl[[p]][[n]] <- cov_list[[n]]
    }else{
      new_cl[[p]][[n]] <- cov_list[[n]][[p]]
    }
  }
  new_cl[[p]] <- brick(new_cl[[p]])
}


## get space-time-locs grid to predict onto
f_orig <- data.table(cbind(coordinates(simple_raster), t=1))
# add time periods
fullsamplespace <- copy(f_orig)
if(nperiods > 1){
  for(p in 2:nperiods){
    tmp <- f_orig
    tmp[,t := p]
    fullsamplespace <- rbind(fullsamplespace,tmp)
  }
}


## get surface to project on to
pcoords <- cbind(x=fullsamplespace$x, y=fullsamplespace$y)
groups_periods <- fullsamplespace$t

## use inla helper functions to project the spatial effect.
A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)

## extract cell values  from covariates, deal with timevarying covariates here
cov_vals <- list()
for(p in 1:nperiods){
  cov_vals[[p]] <- raster::extract(new_cl[[p]], pcoords[1:(nrow(fullsamplespace)/nperiods),])
  cov_vals[[p]] <- (cbind(int = 1, cov_vals[[p]]))
}


#########
#########
## TMB ##
#########
#########


###########
## SETUP ##
###########

if(use_sim_data) X_xp = as.matrix(cbind(1, dt[,covs[, name], with=FALSE]))
if(use_real_data) X_xp = as.matrix(rep(1, nrow(dt)), ncol = 1, nrow = nrow(dt))


Data = list(num_i=nrow(dt),                 ## Total number of observations
            num_s=mesh_s$n,                 ## Number of vertices in SPDE mesh
            num_t=nperiods,                 ## Number of periods
            num_z=1,
            y_i=dt[, Y], ##                 ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            n_i=dt[, N], ##                 ## Number of observed exposures in the cluster (N in binomial likelihood)
            t_i=as.numeric(as.factor(dt[, year]))-1, ## Sample period of ( starting at zero )
            w_i = dt[, weight], 
            X_ij=X_xp,                      ## Covariate design matrix
            M0=spde$param.inla$M0,          ## SPDE sparse matrix
            M1=spde$param.inla$M1,          ## SPDE sparse matrix
            M2=spde$param.inla$M2,          ## SPDE sparse matrix
            Aproj = A.proj,                 ## mesh to prediction point projection matrix
            flag = 1, ##                    ## do normalization outside of optimization if 1
            options = c(1, 1, 0))              ## option1==1 use priors, option2==1 use nugget, option3==1, use adreprt

## staring values for parameters
Parameters = list(alpha_j   =  rep(0,ncol(X_xp)),                 ## FE parameters alphas
                  logtau=1.0,                                     ## log inverse of tau  (Epsilon)
                  logkappa=0.0,	                                  ## Matern Range parameter
                  trho=0.5,
                  zrho=0.5,
                  Epsilon_stz=array(0, dim = c(mesh_s$n, ncol=nperiods)), ## random effects at GP vertex locations
                  nug_i = rep(0, nrow(dt)))     


#########
## FIT ##
#########

templ <- "model"
setwd("/homes/azimmer/tmb_transition/realistic_sims")
TMB::compile(paste0('./', templ,".cpp"))
dyn.load( dynlib(templ) )

## TODO: could also do a simple run to start to get better starting params
## Report0 = obj$report() 

obj <- MakeADFun(data=Data, parameters=Parameters, map=list(zrho = factor(NA)), 
                 random="Epsilon_stz", hessian=TRUE, DLL=templ)
obj <- MakeADFun(data=Data, parameters=Parameters, map=list(zrho = factor(NA),
                                                            trho = factor(NA)),
                 random="Epsilon_stz", hessian=TRUE, DLL=templ)
obj <- normalize(obj, flag="flag")


## Run optimizer
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr))
                        ##    control     =    list(eval.max=1e4, iter.max=1e4, trace=1)))
fit_time_tmb<- proc.time()[3] - ptm

## Get standard errors
SD0 = TMB::sdreport(obj,getJointPrecision=TRUE)
tmb_total_fit_time <- proc.time()[3] - ptm 
tmb_sdreport_time <-  tmb_total_fit_time - fit_time_tmb

#############
## PREDICT ##
#############

message('making predictions')

## now we can take draws and project to space-time raster locs
mu <- c(SD0$par.fixed,SD0$par.random)

## simulate draws
ptm2 <- proc.time()[3]
rmvnorm_prec <- function(mu, prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- Cholesky(prec, super=TRUE)
  z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}
draws <- rmvnorm_prec(mu = mu , prec = SD0$jointPrecision, n.sims = ndraws)
tmb_get_draws_time <- proc.time()[3] - ptm2

## separate out the draws
parnames <- c(names(SD0$par.fixed), names(SD0$par.random))
epsilon_draws <- draws[parnames=='Epsilon_stz',]
alpha_draws   <- draws[parnames=='alpha_j',]

## values of S at each cell (long by nperiods)
cell_s <- as.matrix(A.pred %*% epsilon_draws)

## covariate values by alpha draws
tmb_vals <- list()
for(p in 1:nperiods)  tmb_vals[[p]] <- cov_vals[[p]] %*% alpha_draws

cell_l <- do.call(rbind,tmb_vals)

## add together linear and st components
pred_tmb <- cell_l + cell_s

## save prediction timing
totalpredict_time_tmb <- proc.time()[3] - ptm2

##########
## SAVE ##
##########

## make some useful files first
summ_tmb <- cbind(median = (apply(pred_tmb,1,median)),
                  sd = (apply(pred_tmb,1,sd)))

ras_med_tmb  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,1]),
                                          t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")
ras_sdv_tmb  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_tmb[,2]),
                                          t=rep(1:nperiods,each=nrow(pred_tmb)/nperiods)),"p","t")

saveRDS(file = sprintf('%s/modeling/tmb/outputs/tmb_preds_median_raster.rds', out.dir), object = ras_med_tmb)
saveRDS(file = sprintf('%s/modeling/tmb/outputs/tmb_preds_stdev_raster.rds', out.dir), object = ras_sdv_tmb)


##########
##########
## INLA ##
##########
##########


#########
## FIT ##
#########
A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = dt.coords,
  group = dt.pers
)
space   <- inla.spde.make.index("space",
                                n.spde = spde$n.spde,
                                n.group = nperiods)

inla.covs <- covs$name
design_matrix <- data.frame(int = rep(1, nrow(dt)), dt[, inla.covs, with=F])
stack.obs <- inla.stack(tag='est',
                        data=list(Y=dt$Y), ## response
                        A=list(A,1), ## proj matrix, not sure what the 1 is for
                        effects=list(
                          space,
                          design_matrix)
                        )

formula <- formula(paste0('Y ~ -1+int+',
(paste(inla.covs, collapse = ' + ')),
' + f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))'))

ptm <- proc.time()[3]

inla.setOption("enable.inla.argument.weights", TRUE)
res_fit <- inla(formula,
                data = inla.stack.data(stack.obs),
                control.predictor = list(A = inla.stack.A(stack.obs),
                                         link = 1,
                                         compute = FALSE),
                control.fixed = list(expand.factor.strategy = 'inla'),
                control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                control.compute=list(config = TRUE),
                family = 'binomial',
                num.threads = ncores, ##TODO 
                Ntrials = dt$N,
                weights = rep(1, nrow(dt)),
                verbose = TRUE,
                keep = FALSE)
fit_time_inla <- proc.time()[3] - ptm



#############
## PREDICT ##
#############

ptm <- proc.time()[3]
draws <- inla.posterior.sample(ndraws, res_fit)
inla_get_draws_time <- proc.time()[3] - ptm

## get parameter names
par_names <- rownames(draws[[1]]$latent)

## index to spatial field and linear coefficient samples
s_idx <- grep('^space.*', par_names)
l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor', par_names))

## get samples as matrices
pred_s <- sapply(draws, function (x) x$latent[s_idx])
pred_l <- sapply(draws, function (x) x$latent[l_idx])
rownames(pred_l) <- res_fit$names.fixed


## get samples of s for all coo locations
s <- as.matrix(A.pred %*% pred_s)

## extract cell values  from covariates, deal with timevarying covariates here
inla_vals <- list()
for(p in 1:nperiods)  inla_vals[[p]] <- cov_vals[[p]] %*% pred_l

l <- do.call(rbind,inla_vals)

pred_inla <- s+l

## make them into time bins
len = nrow(pred_inla)/nperiods
totalpredict_time_inla <- proc.time()[3] - ptm

##########
## SAVE ##
##########

## make some useful files first
summ_inla <- cbind(median = (apply(pred_inla,1,median)),
                   sd = (apply(pred_inla,1,sd)))

ras_med_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,1]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")
ras_sdv_inla  <- rasterFromXYZT(data.table(pcoords,p=plogis(summ_inla[,2]), t=rep(1:nperiods,each=nrow(pred_inla)/nperiods)),"p","t")

saveRDS(file = sprintf('%s/modeling/inla/outputs/inla_preds_median_raster.rds', out.dir), object = ras_med_inla)
saveRDS(file = sprintf('%s/modeling/inla/outputs/inla_preds_stdev_raster.rds', out.dir), object = ras_sdv_inla)



################
################
## VALIDATION ##
################
################

dir.create(sprintf('%s/validation', out.dir))

###################################
## summarize fitted param values ##
###################################

res <- data.table(st_mesh_nodes = rep(nrow(epsilon_draws),2))
res[,cores           := rep(ncores,2)]
res[,s_mesh_max_edge := rep(maxedge,2)]
res[,periods         := c(4,4)]
res[,draws           := c(ndraws,ndraws)]

## time variables
res[,fit_time  := c(fit_time_inla,fit_time_tmb)]
res[,pred_time := c(totalpredict_time_inla,totalpredict_time_tmb)]
res[,pt_tmb_sdreport_time := c(NA,tmb_sdreport_time)]
res[,pt_get_draws_time := c(inla_get_draws_time,tmb_get_draws_time)]

## fe coefficients
for(i in 1:length(res_fit$names.fixed)){
  fn <- res_fit$names.fixed[i]
  res[[paste0('fixedeff_',fn,'_mean')]] <- c(res_fit$summary.fixed$mean[i], SD0$par.fixed[i])
  res[[paste0('fixedeff_',fn,'_sd')]]   <- c(res_fit$summary.fixed$sd[i], sqrt(SD0$cov.fixed[i,i]))
}

## hyperparameters
res[,hyperpar_logtau_mean := c(res_fit$summary.hyperpar[1,1], SD0$par.fixed['logtau']) ]
res[,hyperpar_logtau_sd := c(res_fit$summary.hyperpar[1,2], sqrt(SD0$cov.fixed['logtau','logtau'])) ]

res[,hyperpar_logkappa_mean := c(res_fit$summary.hyperpar[2,1], SD0$par.fixed['logkappa']) ]
res[,hyperpar_logkappa_sd := c(res_fit$summary.hyperpar[2,2], sqrt(SD0$cov.fixed['logkappa','logkappa'])) ]

res[,hyperpar_rho_mean := c(res_fit$summary.hyperpar[3,1], SD0$par.fixed['trho']) ]
res[,hyperpar_rho_sd := c(res_fit$summary.hyperpar[3,2], sqrt(SD0$cov.fixed['trho','trho'])) ]

## truth
true.params <- c(rep(NA, 9),
                 alpha, NA, ## intercept
                 c(rbind(betas, rep(NA, length(betas)))), ## fixed effect coefs
                 -0.5 * log(4 * pi * sp.var * sp.kappa^2), NA, ## log tau
                 log(sp.kappa), NA, ## log kappa
                 t.rho, NA ## time rho
                 )

rr <- data.table(item=colnames(res))
rr <- cbind(rr,true.params, t(res))
names(rr) <- c('quantity','TRUE', 'R-INLA','TMB')
rr$diff <- rr[,3]-rr[,4]

## we can now plot this table with: grid.table(rr)

####################
## setup big plot ##
####################
pdf(sprintf('%s/validation/inla_tmb_summary_comparison_plot.pdf',out.dir), height=15,width=30)

grid.table(rr)

plot.in.logit.space <- FALSE

for(thing in c('median','stdev')){
    layout(matrix(1:20, 4, 5, byrow = TRUE))
    samp=sample(cellIdx(ras_med_inla[[1]]),1e4)
    for(i in 1:4){

      ## TODO: allow plotting in logit space
      if(plot.in.logit.space){
        
      }else{
        
        if(thing=='median'){
          rinla <-  ras_med_inla[[i]]
          rtmb  <-  ras_med_tmb[[i]]
        }
        if(thing=='stdev'){
          rinla <-  ras_sdv_inla[[i]]
          rtmb  <-  ras_sdv_tmb[[i]]
        }

##        true <- true.rast[[i]] ## in logit space
##        values(true) <- plogis(values(true))
      }

      
      tmp <- subset(dt, period_id==i)
      
      # rasters
      par(mar = c(0, 0, 1.4, 1),bty='n')
      maxes <- max(c(as.vector(rtmb),as.vector(rinla)),na.rm=TRUE)
      plot(rinla-rtmb, maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend=T, main=paste0('DIFFERENCE ',thing))
      plot(simple_polygon, add = T)
      plot(rtmb,  maxpixel=1e7, col=rainbow(100), axes=FALSE, legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main='TMB',zlim=c(0,maxes))
      plot(simple_polygon, add = T)
      plot(rinla, maxpixel=1e7, col=rainbow(100), axes=FALSE, legend=FALSE, main='R-INLA',zlim=c(0,maxes))
      plot(simple_polygon, add = T)
##      plot(true, maxpixel=1e7, col=rainbow(100), axes=FALSE, legend=FALSE, main='TRUE (mean)',zlim=c(0,maxes))
      
      # scatter
      par(mar = c(4, 4, 2, 2),bty='n')
      plot(x=as.vector(rinla)[samp],y=as.vector(rtmb)[samp],xlab='R-INLA',ylab='TMB',cex=.01,pch=19,main='COMPARE')
      lines(x=c(0,maxes),y=c(0,maxes),col='red')

      # plot data points over the shapefile
      plot(simple_polygon, main='DATA LOCATIONS')
      points( x=tmp$long,y=tmp$lat, pch=19, cex=(tmp$N / max(tmp$N)) )

      # residual
      #  tmp$inla<-extract(ras_med_inla[[i]],cbind(tmp$longitude,y=tmp$latitude))
      #  tmp$tmb<-extract(ras_med_tmb[[i]],cbind(tmp$longitude,y=tmp$latitude))
      #  tmp$dat <- tmp$died/tmp$N
      #  tmp$resid_inla <- tmp$dat-tmp$inla
      #  tmp$resid_tmb  <- tmp$dat-tmp$tmb
      #  tmp<-subset(tmp,dat<quantile(tmp$dat,.9))
      #  plot(x=tmp$dat,y=tmp$resid_inla, pch=19,col='red',cex=.1,main='RESIDUALS')
      #  points(x=tmp$dat,y=tmp$resid_tmb, pch=19,col='blue',cex=.1)
    }
}

layout(matrix(1, 1, 1, byrow = TRUE))
####
#### Compare mean and distribution of random effects
summ_gp_tmb  <- t(cbind((apply(epsilon_draws,1,quantile,probs=c(.1,.5,.9)))))
summ_gp_inla <- t(cbind((apply(pred_s,1,quantile,probs=c(.1,.5,.9)))))
  # all time-space random effects

plot_d <- data.table(tmb_median = summ_gp_tmb[,2],inla_median = summ_gp_inla[,2],
                     tmb_low    = summ_gp_tmb[,1],inla_low    = summ_gp_inla[,1],
                     tmb_up     = summ_gp_tmb[,3],inla_up     = summ_gp_inla[,3])

plot_d$period <- factor(rep(1:4,each=nrow(plot_d)/4))
plot_d$loc    <- rep(1:(nrow(plot_d)/4),rep=4)

if(nrow(plot_d)>2500)
  plot_d <- plot_d[sample(nrow(plot_d),2500,replace=F),]


ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
  geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
  ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

# plot locations where they are different, are they near or far from data?
plot_d[, absdiff := abs(tmb_median-inla_median)]
nodelocs <- do.call("rbind", replicate(4, mesh_s$loc, simplify = FALSE))
biggdiff <- unique(nodelocs[which(plot_d$absdiff>quantile(plot_d$absdiff,prob=0.80)),])

nodelocs <- cbind(nodelocs,plot_d)
if(nrow(nodelocs)>2500)
  nodelocs <- nodelocs[sample(nrow(nodelocs),2500,replace=FALSE),]

par(mfrow=c(2,2))
for(i in 1:4){
  plot(simple_polygon, main='Mesh nodes sized by abs difference TMB and R-INLA')
  points(x=dt$longitude[dt$period_id==i],y=dt$latitude[dt$period_id==i], pch=19, cex=0.1)
  points(x=nodelocs$V1[nodelocs$period==i],y=nodelocs$V2[nodelocs$period==i], pch=1, cex=nodelocs$absdiff[nodelocs$period==i]*5, col='red')
}

# catterpillar plot
plot_d <- plot_d[order(period,tmb_median)]
plot_d[,i := seq(1,.N), by = period]
ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
  geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='blue', size=.8, alpha=.3) +
  geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='red', size=.8, alpha=.3) +
  facet_wrap(~period) +
  ggtitle('Comparison of random effects (10% to 90% quantiles) ... RED == R-INLA ... BLUE == TMB')


dev.off()


    


