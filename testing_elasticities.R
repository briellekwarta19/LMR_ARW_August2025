#install.packages("devtools")
#devtools::install_github("AldridgeCaleb/meta-pop-bio")
library("metapopbio")
library(tidyverse)
#loading data
path <- here::here()
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))

c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

#here I am making a matrix where each value maps to the correct patch/stage/rate 
patch_names <- rep(1:n_patches, each = n_stages)
stage_names <- rep(c('juvenile', 'adult'), times = n_patches)


mat.names <- matrix(NA, nrow = n_patches*n_stages, ncol = n_patches*n_stages)
rownames(mat.names) <- patch_names
colnames(mat.names) <- stage_names
mat.names[seq(1,nrow(mat.names),by = 2),] <- 'fecund'
mat.names[seq(2,nrow(mat.names),by = 2),] <- 'surv'

mat.id <- matrix(NA, nrow = n_patches*n_stages, ncol = n_patches*n_stages)
for(i in 1:nrow(mat.id)){
  for(j in 1:ncol(mat.id)){
    mat.id[i,j] <- paste0('patch: ', rownames(mat.names)[i], ', rate: ', colnames(mat.names)[j], ' ', mat.names[i,j] )
  }
}

#test [1,1]
mat.id[1:2,1:2]

P <- metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by)

#-------------------------------------------------------------------------------
#Elasticity analysis: 

#1. Lets do this first for the basemodel
A_base <- spmm.project.matrix(P, #vec-permutation matrix,
                         BB, #block diagonal matrix for demographic paramters 
                         MM, #block diagonal for movement paramters
                         group_by = group_by, #grouping projections (by patches here)
                         lh_order = lh_order) #order of events (demographic before movement here)


#extracting elasticities for the base model (no harvest, no deterrents)
BB_e_base <- spmm.demo.elas(BB,A_base,P,MM)

#extracting location of top 5 rates
top5_BB_base <- order(BB_e_base, decreasing = TRUE)[1:5]
#extracting indices of those rates
top5_BB_base_indices <- arrayInd(top5_BB_base, .dim = dim(BB_e_base))

#now combine into a dataframe which matches with the rate name:
top5_BB_base_df <- data.frame(rate = mat.id[top5_BB_base_indices], value = BB_e_base[top5_BB_base])
top5_BB_base_df$type <- 'base model'


#2. Lets do this first for a different strategy: e.g., harvest 0.1 at first 10 patches
harv_mat <- unblk.diag(BB, n_stages)

#Need to make a new BB matrix for the new harvest
for (i in 1:10) {
  B <- harv_mat[i]
  
  M <- -log(B[[1]][-1, ])  
  M <- M + 0.5 #This will be adjusted per strategy
  
  B[[1]][-1, ] <- exp(-M)
  harv_mat[i] <- B
}

BB_harv <- blk.diag(harv_mat)

#check that BB_harv is different from BB
all.equal(BB_harv, BB)

A_harv <- spmm.project.matrix(P, #vec-permutation matrix,
                              BB = BB_harv, #block diagonal matrix for demographic paramters 
                              MM, #block diagonal for movement paramters
                              group_by = group_by, #grouping projections (by patches here)
                              lh_order = lh_order) #order of events (demographic before movement here)


#extracting elasticities for this strategy
BB_e_harv <- spmm.demo.elas(BB_harv,A_harv,P,MM)

#find elasticity values from the harvest BB model that match the top 5 ones in the base model
#and combine into a dataframe which matches with the rate name:
top5_BB_harv_df <- data.frame(rate = mat.id[top5_BB_base_indices], value = BB_e_harv[top5_BB_base])
top5_BB_harv_df$type <- 'harvest strategy'


#3. combine top5 dataframes
top5_bb <- rbind(top5_BB_base_df, top5_BB_harv_df)

ggplot(top5_bb, aes(x = value, y = rate, fill = type))+
  geom_col(position = "dodge") + ggtitle("BB elasticities")

#Question: should BB_e_harv equal BB_e_base??
#all.equal(BB_e_harv, BB_e_base)





