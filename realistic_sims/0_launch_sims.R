## this script can be used to launch 1_run_simulation.R in parallel on the IHME cluster
## written by aoz


####################
## run some setup ##
####################

################################################################################
logging_note <- 'parallel test run'
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

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

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

###############################
## setup things to loop over ##
###############################

## this list should align with the args in 1_run_simulation.R

reg <- 'nga'## commandArgs 4
year_list <- 2000 ## commandArgs 5
cov_names <- "c('access2', 'distrivers', 'evi'   , 'mapincidence')" 
cov_measures <- "c('mean'   , 'mean'      , 'median', 'mean')"
betas <- "c(.5, -1, 1, -.5)"
alpha <- 0
sp.range <-  sqrt(8)     ## commandArgs 10 sqrt(8)  ## kappa=sqrt(8)/sp.range, so sp.range=sqrt(8) -> kappa=1 -> log(kappa)=0 (for R^2 domain)
sp.var <- 0.5            ## sp.var = 1/(4*pi*kappa^2*tau^2) (for R^2 domain)
sp.alpha <- 2.0          ## matern smoothness = sp.alpha - 1 (for R^2 domain)
nug.var <- .5 ^ 2        ## nugget variance
t.rho <-  0.8            ## annual temporal auto-corr
mesh_s_max_edge <- c("c(0.2,5)",
                     "c(0.3,5)") ## commandArgs 15
n.clust <-  50                   ## clusters PER TIME slice
m.clust <- 35                    ## mean number of obs per cluster (poisson)
sample.strat <- 'random'         ## random or by population for now. ## TODO cluser design
cores <- 5
ndraws <- 250                    ## commandArgs 20

loopvars <- expand.grid(reg,
                        year_list,
                        cov_names,
                        cov_measures,
                        betas,
                        alpha,
                        sp.range,
                        sp.var,
                        sp.alpha,
                        nug.var,
                        t.rho,
                        mesh_s_max_edge,
                        n.clust,
                        m.clust,
                        sample.strat,
                        cores,
                        ndraws)

loopvars$run_date <- NA ## keep track of run_dates to later compare runs

## prepare a set of run_dates so we can write the complete loopvars to each run_date dir
for(ii in 1:nrow(loopvars)){
  loopvars$run_date[ii] <- make_time_stamp(TRUE)
  Sys.sleep(1.1)
}

for(ii in 1:nrow(loopvars)){

  ## make a run_date and setup output directory
  run_date <- loopvars$run_date[ii]

  ## now we can setup our main directory to save these results and log our note and stdouts
  out.dir  <- sprintf('/homes/azimmer/tmb_inla_sim/%s', run_date)
  dir.create(out.dir)
  dir.create(paste0(out.dir, '/errors'))
  dir.create(paste0(out.dir, '/output'))
  
  ## write the log note
  fileConn <- file(sprintf("%s/run_notes.txt", out.dir))
  writeLines(logging_note, fileConn)
  close(fileConn)

  ## save loopvars to this dir to reload into the parallel env
  saveRDS(file = paste0(out.dir, '/loopvars.rds'), obj = loopvars)

  ## save and reload loopvars in parallel env. that way, we only need to pass in iter/row #
  qsub.string <- qsub_sim(iter = ii, ## sets which loopvar to use in parallel
                          run_date = run_date,
                          slots = 4, 
                          codepath = '/homes/azimmer/tmb_transition/realistic_sims/1_run_simulation.R', 
                          singularity = 'default',
                          singularity_opts = NULL,
                          logloc = NULL ## defaults to input/output dir in sim run_date dir
                          )

  ## launch the job
  
  system(qsub.string)
}
