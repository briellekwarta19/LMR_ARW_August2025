# Set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load packages -----------------------------------------------------------
if (!require("pacman")) {
  install.packages(pacman)
  library("pacman")
}
p_load(metapopbio)


# Load data ---------------------------------------------------------------
path <- getwd()
filename <- "/LMRARW-ICPMM-(SPMM).xlsm"
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
P <- metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by)
A <- spmm.project.matrix(P, BB, MM, 
                         group_by = group_by, 
                         lh_order = lh_order)
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
spmm.plot(projs, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = patch_names)


# Construct utility -------------------------------------------------------
## Function to calculate cost
calculate.cost <- function(harvest_mortality) {
  initial_cost <- 100000  # Initial cost for 0.01 harvest mortality
  growth_rate <- 1.1  # 10% increase
  cost <- initial_cost * (growth_rate) ^ (harvest_mortality / 0.01)
  return(cost)
}

## Construct function to calculate utility
final.p <- function(xx) {
  x1 <- xx[1]
  projs <-
    spmm.project(
      n = n,
      A = A,
      BB = BB,
      MM = MM,
      P = P,
      n_timesteps = n_timesteps,
      n_stages = n_stages,
      n_patches = n_patches,
      harv = c(x1)
    )
  final_p <- (1 - ((sum(projs[, n_timesteps]) / n_patches) / 10000))
  cost <- calculate.cost(x1)
  final_cost <- (1 - (cost / 20000000))
  
  ## Equal objective weights
  U <- final_p * 0.5 + final_cost * 0.5
  return(U)
}


# Optimize harvest strategy using genetic algorithm search ----------------
p_load(rgenoud)

F_est <- genoud(
  final.p,
  nvars = 1,
  max = TRUE,
  Domains = matrix(c(0.01, 0.2), ncol = 2, byrow = TRUE),
  boundary.enforcement = 2
)

## Projection using optimal harvest level
projs <-
  spmm.project(
    n = n,
    A = A,
    BB = BB,
    MM = MM,
    P = P,
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches,
    harv = F_est$par
  )
spmm.plot(projs, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = patch_names)


# Notes -------------------------------------------------------------------
# Need to add in options for deterrent (see deter arg in ?spmm.project)

# Could model the random proportion of the population that moves rbinom(1, 1, p)

# When estimating rel. ab. / biomass should sampling be stratified? 
# E.g., LTRM: Main channel, side channel, backwater, impounded shoreline
# Others: embayment, oxbow, tributary

## Function for proportional scoring
prop.scaling <- function(x = NULL, direction = "max") {
  low <- min(x, na.rm = TRUE)
  high <- max(x, na.rm = TRUE)
  if (direction == "max") {
    y <- (x - low) / abs(low - high)
    return(y)
  } if (direction == "min") {
    y <- (high - x) / abs(low - high)
    return(y)
  }
}
