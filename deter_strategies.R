#install.packages("devtools")
#devtools::install_github("AldridgeCaleb/meta-pop-bio")
library("metapopbio")

path <- here::here()
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))
## Assign objects from list
c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

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


