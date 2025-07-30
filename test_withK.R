
path <- here::here()
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(demo <- spmm.readxl(path, filename))

filename2 <- here::here("nid-info.xlsx")
nid_info <- read.xlsx(filename2)


## Assign objects from list
c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% demo

## "Manual" projection and plotting
P <- metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by)

A <- spmm.project.matrix(P, BB, MM, 
                         group_by = group_by, 
                         lh_order = lh_order)

K <- nid_info$K

projs <- spmm.project(
  n = n,  # number of stage/age animals in patch i
  A = A,
  BB = BB,
  MM = MM,
  P = P,
  ddf = list(K = K, beta = 0.2, theta = 1.1),
  n_timesteps = n_timesteps,
  n_stages = n_stages,
  n_patches = n_patches, 
  mod_mort = 0.2
)

spmm.plot(projs, 
          ylabs = "Biomass",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = patch_names,
          ylim_max = "patch")
