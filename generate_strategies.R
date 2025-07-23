library("metapopbio")
library(plyr)

#Load data ---------------------------------------------------------------
path <- here::here()
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))
## Assign objects from list
c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

#Create data frame that lists strategy names:
harv_vals <- seq(0,0.2,0.05)
deter_names <-  c('None', 'One lower', 'Two lower', 'One leading edge',
                  'Two leading edge', 'One lower + one leading edge ', 'Two lower + one leading edge',
                  'Two lower + two leading edge')

deter_id <- seq(1:8)

strategy_names <- expand.grid(harv = harv_vals, deter = deter_id)


#Create the vec-permutation matrix (dimension = stages*patches x stages*patches)
#this stays the same across all management strategies
P <- metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by)

#different MM levels:
base_movement <- unblk.diag(MM, n_patches)

modify.movement <- function(movement_matrix, position, lower_factor = 0.95, upper_factor = 0.5) {
  lwr_idx <- seq_len(position - 1)
  upr_idx <- if (position < nrow(movement_matrix)) (position + 1):nrow(movement_matrix) else integer(0)
  movement_matrix[lwr_idx, position] <- movement_matrix[lwr_idx, position] * lower_factor
  movement_matrix[upr_idx, position] <- movement_matrix[upr_idx, position] * upper_factor
  movement_matrix[position, position] <- NA
  movement_matrix[position, position] <- 1 - sum(movement_matrix[-position, position], na.rm = TRUE)
  movement_matrix[, position] <- movement_matrix[, position] / sum(movement_matrix[, position])
  return(movement_matrix)
}


### 2. One lower ("cut off spigot")
move_mod_2 <- base_movement
pos <- 2
for (i in 1:n_stages) {
  move_mod_2[[i]] <- modify.movement(move_mod_2[[i]], pos)
}
MM_2 <- blk.diag(move_mod_2)

### 3. Two lower ("cut off spigot and main")
move_mod_3 <- base_movement
pos <- 2:3
for (i in 1:n_stages) {
  for (j in pos) {
    move_mod_3[[i]] <- modify.movement(move_mod_3[[i]], j)
  }
}
MM_3 <- blk.diag(move_mod_3)

### 4. One leading edge ("line in the sand")
move_mod_4 <- base_movement
pos <- 9
for (i in 1:n_stages) {
  move_mod_4[[i]] <- modify.movement(move_mod_4[[i]], pos)
}
MM_4 <- blk.diag(move_mod_4)

### 5. Two leading edge ("two lines in the sand")
move_mod_5 <- base_movement
pos <- 9:10
for (i in 1:n_stages) {
  for (j in pos) {
    move_mod_5[[i]] <- modify.movement(move_mod_5[[i]], j) 
  }
}
MM_5 <- blk.diag(move_mod_5)

### 6. One lower, one leading edge ("squeeze / sandwich") 
move_mod_6 <- base_movement
pos <- c(2, 9)
for (i in 1:n_stages) {
  for (j in pos) {
    move_mod_6[[i]] <- modify.movement(move_mod_6[[i]], j)
  }
}
MM_6 <- blk.diag(move_mod_6)

### 7. Two lower, one leading edge ("mcdouble")
move_mod_7 <- base_movement
pos <- c(2:3, 9)
for (i in 1:n_stages) {
  for (j in pos) {
    move_mod_7[[i]] <- modify.movement(move_mod_7[[i]], j)
  }
}
MM_7 <- blk.diag(move_mod_7)

### 8. Two lower, two leading edge ("big mac")
move_mod_8 <- base_movement
pos <- c(2:3, 9:10)
for (i in 1:n_stages) {
  for (j in pos) {
    move_mod_8[[i]] <- modify.movement(move_mod_8[[i]], j)
  }
}
MM_8 <- blk.diag(move_mod_8)

#-------------------------------------------------------------------------------
#### Simulate strategies ####

MMs <- list(MM, MM_2, MM_3, MM_4, MM_5, MM_6, MM_7, MM_8)

harv.strategies <- array(0, c(length(strategy_names$harv), n_patches))
for(i in 1:length(strategy_names$harv)){
  harv.strategies[i,1:10] <- strategy_names$harv[i]
}

A <- list()
projs <- list()

for(i in 1:length(strategy_names$harv)){
  A[[i]] <- spmm.project.matrix(P, 
                           BB, 
                           MMs[[strategy_names$deter[i]]], 
                           group_by = group_by, #change to grouping
                           lh_order = lh_order)
   

  projs[[i]] <- spmm.project(
    n = n,  # number of stage/age animals in patch i
    A = A[[i]],
    BB = BB,
    MM = MMs[[strategy_names$deter[i]]],
    P = P,
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches, 
    mod_mort = harv.strategies[i,]
  )
}


#### calculate cost per strategy ####
calculate.cost <- function(harvest_mortality) {
  initial_cost <- 100000  # Initial cost for 0.01 harvest mortality
  growth_rate <- 1.1  # 10% increase
  cost <- ifelse(harvest_mortality == 0, 0, initial_cost * (growth_rate) ^ (harvest_mortality / 0.01))
  return(cost)
}

# deter_names <-  c('None', 'One lower', 'Two lower', 'One leading edge',
#                   'Two leading edge', 'One lower + one leading edge ', 'Two lower + one leading edge',
#                   'Two lower + two leading edge')

deter_nums <-  c(0, 1, 2, 1, 2, 2, 3,4) #number of deterrants per deterrant strategy

strategy_cost <- expand.grid(harv = harv_vals, deter = deter_nums)


strategy_cost <- strategy_cost %>% mutate(
  cost = calculate.cost(harv) + deter*calculate.cost(0.2)
)

### checking out outcomes of objectives:
final_pop <- final_dist <- rep(NA, length(projs))

for(i in 1:length(projs)){
  projs_select <- projs[[i]]
  final_pop[i] <- sum(projs_select[, n_timesteps])
  
  summed_pop <- projs_select[seq(1, nrow(projs_select), by = 2), ] + projs_select[seq(2, nrow(projs_select), by = 2), ]

  final_dist[i] <- sum(summed_pop[,n_timesteps] >= 0)
}

strategy_outcomes <- strategy_cost

strategy_outcomes$deter <- deter_names[strategy_names$deter]
strategy_outcomes$pop <- final_pop #final population
strategy_outcomes$num.invaded.patches <- final_dist # of patches with population

#checking out outcomes for top harvest strategy:

strategy_outcomes %>% arrange(pop)
strategy_outcomes %>% filter(harv == 0.2) %>% arrange(pop)
