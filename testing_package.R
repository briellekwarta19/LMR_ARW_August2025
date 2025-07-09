# Set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load packages -----------------------------------------------------------
if (!require("pacman")) {
  install.packages(pacman)
  library("pacman")
}
p_load(metapopbio)


# Load data ---------------------------------------------------------------
path <- here::here("LMRARW-ICPMM-AR-Demo-Files")
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))
## Assign objects from list
c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat


# Project and plot --------------------------------------------------------
## Automated function
auto.spmm(
  path = path,
  filename = filename,
  plot = TRUE
)

## "Manual" projection and plotting
#1. creating the vec-permutation matrix (dimension = stages*patches x stages*patches)
P <- metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by)

#2. Create projection matrix
A <- spmm.project.matrix(P, #vec-permutation matrix,
                         BB, #block diagonal matrix for demographic paramters 
                         MM, #block diagonal for movement paramters
                         group_by = group_by, #grouping projections (by patches here)
                         lh_order = lh_order) #order of events (demographic before movement here)

#3. Make projections
projs <- spmm.project(
  n = n,  # number of stage/age animals in patch i
  A = A,
  BB = BB,
  MM = MM,
  P = P,
  n_timesteps = n_timesteps,
  n_stages = n_stages,
  n_patches = n_patches
)

#plotting everything:
spmm.plot(projs, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = patch_names)


#testing to plot just a few patches:
subset <- projs[1:2,]
comment(subset) <- comment(projs)

spmm.plot(subset, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = patch_names[1])

#making it generalizable
patch_match <- rep(patch_names, each = n_stages)

ex_patch <- 'MS River & Lower L&Ds'

subset <- projs[which(patch_match == ex_patch),]
comment(subset) <- comment(projs)

spmm.plot(subset, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = ex_patch)

#### adding in Harvest ####
projs_harvest <- spmm.project(
  n = n,  # number of stage/age animals in patch i
  A = A,
  BB = BB,
  MM = MM,
  P = P,
  n_timesteps = n_timesteps,
  n_stages = n_stages,
  n_patches = n_patches, 
  harv = 0.2
  #harv = rep(0.2, nrow(BB))
  #harv = rep(0.2, nrow(BB)* ncol(BB))
  #harv = matrix(0.2, nrow = nrow(BB), ncol =  ncol(BB))
  #harv = rep(0.2, ((nrow(BB)* ncol(BB)) -ncol(BB)) / 2)
  #harv = rep(0.2, (nrow(BB)* ncol(BB)) -ncol(BB))
)

ex_patch <- 'MS River & Lower L&Ds'

subset_harv <- projs[which(patch_match == ex_patch),]
comment(subset_harv) <- comment(projs_harvest)

spmm.plot(subset_harv, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = ex_patch)

spmm.plot2(subset_harv, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = ex_patch)



