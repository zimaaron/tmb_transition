## this script can be used to launch 1_run_simulation.R in parallel on the IHME cluster
## written by aoz


####################
## run some setup ##
####################

################################################################################
logging_note <- 'Running these models with real HIV SSSA data. using weights'
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
year_list <- 2000, ## commandArgs 5
cov_names <- "c()",
cov_measures,
betas,
alpha,
sp.range, ## commandArgs 10
sp.alpha,
nug.var,
t.rho,
mesh_s_max_edge,
n.clust, ## commandArgs 15
m.clust,
sample.strat,
cores, 
ndraws,
run_date, ## commandArgs 20
sep = " 
