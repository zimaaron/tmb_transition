## this script simulates some realistic datasets for comparison between INLA and TMB
## it leverages existing architecture that the LBD team at IHME has already created
## written by AOZ
## last editted 5/14/18

options(error = recover)

## DO THIS!
################################################################################
## ADD A NOTE! to help identify what you were doing with this run
logging_note <- 'This is a fresh test after turning sims into a function.'
################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## now we can setup our main directory to save these results and log our note
run_date <- make_time_stamp(TRUE)
out.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', run_date)
dir.create(out.dir)
fileConn <- file(sprintf("%s/run_notes.txt", out.dir))
writeLines(logging_note, fileConn)
close(fileConn)

## Now we can switch to the TMB repo
setwd(tmb_repo)
if(pull_tmb_git) system(sprintf('cd %s\ngit pull %s %s', core_repo, remote, branch))
source('./realistic_sims/realistic_sim_utils.R')

## setup is now done. setup some parameters for this simulation

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##############################
## setup tunable parameters ##
##############################

## TODO, make this one big function with these params as args

## pick a country (or region) using iso3 codes (or region names)
reg       <- 'NGA'
year_list <- seq(2000, 2015, by = 5)

## fixed effects betas
betas <- c(.5, -1, 1, -.5)

## covariates
## load some covariates (need matching name and measure)
covs <- data.table(name = c('access2', 'distrivers', 'evi'   , 'mapincidence'),
                   meas = c('mean'   , 'mean'      , 'median', 'mean'))

## gp and ar1 param
sp.kappa <- 1.0
sp.var   <- 0.5
sp.alpha <- 2.0
t.rho    <- 0.8

## number of clusters, mean and sd of cluster size
n.clust <- 100 ## clusters PER TIME slice
m.clust <- 35  ## mean number of obs per cluster (poisson)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################
## load in region/counry shapes and covs ##
###########################################

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
sim.obj <- sim.realistic.data(reg = reg,
                              year_list = year_list,
                              betas = betas,
                              sp.kappa = sp.kappa,
                              sp.alpha = sp.alpha,
                              t.rho = t.rho,
                              n.clust = n.clust,
                              m.clust = m.clust,
                              covs = covs,
                              simple_raster = simple_raster,
                              simple_polygon = simple_polygon, 
                              out.dir = out.dir,
                              seed = 123456)

dt <- sim.obj$sim.dat ## simulated data, lat-long, year, covs, true surface
covs.gp <- sim.obj$cov_gp_raster   ## rasters of covs and true simulated gp field

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################
#####################################
## setup for tmb and INLA modeling ##
#####################################
#####################################
dir.create(sprintf('%s/modeling_inputs/', out.dir))

###############################
## SETUP SOME SHARED OBJECTS ##
###############################

## set some stuff up
dt[, id := 1:.N]
dt.coords <- cbind(dt$long,dt$lat)
nperiod   <- length(unique(dt$year))
dt.pers   <- as.numeric(as.factor(dt[, year]))

## Build spatial mesh over modeling area
## same mesh across time (TODO: only for now)
mesh_s <- inla.mesh.2d(,
                       inla.sp2segment(simple_polygon)$loc,
                       max.edge=c(1,5),
                       cutoff=0.2)
pdf(sprintf('%s/modeling_inputs/mesh.pdf', out.dir))
plot(mesh_s)
plot(covs.gp[[1]], add = TRUE) ## just to show loc of simple_raster under mesh for scale
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
saveRDS(file = sprintf('%s/modeling_inputs/mesh.rds', out.dir), mesh_s)
saveRDS(file = sprintf('%s/modeling_inputs/spde.rds', out.dir), spde)

#########
#########
## TMB ##
#########
#########

###########
## SETUP ##
###########

X_xp = as.matrix(cbind(1, dt[,covs[, name], with=FALSE]))


Data = list(num_i=nrow(dt),                 ## Total number of observations
            num_s=mesh_s$n,                 ## Number of vertices in SPDE mesh
            num_t=nperiod,                  ## Number of periods
            num_z=1,
            y_i=dt[, Y], ##                 ## Number of observed deaths in the cluster (N+ in binomial likelihood)
            n_i=dt[, N], ##                 ## Number of observed exposures in the cluster (N in binomial likelihood)
            t_i=as.numeric(as.factor(dt[, year]))-1, ## Sample period of ( starting at zero )
            w_i=rep(1,nrow(dt)),
            X_ij=X_xp,                      ## Covariate design matrix
            M0=spde$param.inla$M0,          ## SPDE sparse matrix
            M1=spde$param.inla$M1,          ## SPDE sparse matrix
            M2=spde$param.inla$M2,          ## SPDE sparse matrix
            Aproj = A.proj,                 ## mesh to prediction point projection matrix
            options = c(1))                 ## option1==1 use priors


## staring values for parameters
Parameters = list(alpha_j   =  rep(0,ncol(X_xp)),                 ## FE parameters alphas
                  logtau=1.0,                                     ## log inverse of tau  (Epsilon)
                  logkappa=0.0,	                                  ## Matern Range parameter
                  trho=0.5,
                  zrho=0.5,
                  Epsilon_stz=matrix(1, nrow=mesh_s$n, ncol=nperiod))     ## GP locations


#########
## FIT ##
#########

templ <- "model"
TMB::compile(paste0(templ,".cpp"))
dyn.load( dynlib(templ) )

#openmp(10)
obj <- MakeADFun(data=Data, parameters=Parameters,  map=list(zrho=factor(NA)), random="Epsilon_stz", hessian=TRUE, DLL=templ)


## Run optimizer
ptm <- proc.time()[3]
opt0 <- do.call("nlminb",list(start       =    obj$par,
                              objective   =    obj$fn,
                              gradient    =    obj$gr))
                        #      control     =    list(eval.max=1e4, iter.max=1e4, trace=1)))
tmb_fit_time <- proc.time()[3] - ptm


# Get standard errors
## Report0 = obj$report()
ptm <- proc.time()[3]
SD0 <- sdreport(obj,getReportCovariance=TRUE)
## fe_var_covar <- SD0$cov.fixed
tmb_sdreport_time <- proc.time()[3] - ptm


logtau  = SD0$par.fixed['logtau']
logkappa = SD0$par.fixed['logkappa']
unname(sqrt(8.0) / exp(logkappa))
unname(1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa)))



##### Prediction
message('making predictions')
#mu    <- c(SD0$par.fixed[names(SD0$par.fixed)=='alpha'],SD0$par.random[names(SD0$par.random)=="epsilon"])
mu    <- c(SD0$value)

sigma <- SD0$cov

### simulate draws
require(MASS)
npar   <- length(mu)
ndraws <- 50



#############
## PREDICT ##
#############


##########
##########
## INLA ##
##########
##########


#########
## FIT ##
#########


#############
## PREDICT ##
#############


################
################
## VALIDATION ##
################
################
