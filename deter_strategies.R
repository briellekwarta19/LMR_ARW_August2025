library("metapopbio")

path <- here::here()
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))

c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

n_patches <- 14

base_movement <- matrix(c(
  NA
), nrow = n_patches, 
ncol = n_patches, byrow = TRUE)
base_movement[1, ] <- c(0.95, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
base_movement[2, ] <- c(0.05, 0.925, 0.025, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
base_movement[3, ] <- c(0.05, 0.15, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
base_movement[4, ] <- c(0, 0.1, 0.05, 0.85, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
base_movement[5, ] <- c(0, 0, 0, 0, 0.975, 0.025, 0, 0, 0, 0, 0, 0, 0, 0)
base_movement[6, ] <- c(0, 0, 0, 0, 0.01, 0.98, 0.01, 0, 0, 0, 0, 0, 0, 0)
base_movement[7, ] <- c(0, 0, 0, 0, 0, 0.05, 0.95, 0, 0, 0, 0, 0, 0, 0)
base_movement[8, ] <- c(0, 0, 0, 0, 0, 0, 0.05, 0.925, 0.025, 0, 0, 0, 0, 0)
base_movement[9, ] <- c(0, 0, 0, 0, 0, 0, 0, 0.01, 0.99, 0, 0, 0, 0, 0)
base_movement[10, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0.025, 0.95, 0.025, 0, 0, 0)
base_movement[11, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99, 0.01, 0, 0)
base_movement[12, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.98, 0.01, 0)
base_movement[13, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.9, 0)
base_movement[14, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, 0.01, 0.79)

modify.movement <- function(movement_matrix, position, lower_factor = 0.95, upper_factor = 0.5) {
  lwr_idx <- seq_len(position - 1)
  upr_idx <- ifelse(position < ncol(movement_matrix), (position + 1):ncol(movement_matrix), integer(0))
  movement_matrix[position, lwr_idx] <- movement_matrix[position, lwr_idx] * lower_factor
  movement_matrix[position, upr_idx] <- movement_matrix[position, upr_idx] * upper_factor
  movement_matrix[position, position] <- NA
  movement_matrix[position, ] <- ifelse(is.na(movement_matrix[position, ]), 
                                        1 - sum(movement_matrix[position, ], na.rm = TRUE), 
                                        movement_matrix[position, ])
  return(movement_matrix)
}

### 2. One lower ("cut off spigot")
move_mod_2 <- movement
pos <- 2
move_mod_2 <- modify.movement(move_mod_2, pos)

### 3. Two lower ("cut off spigot and main")
move_mod_3 <- movement
pos <- 2:3
for (i in pos) {
  move_mod_3 <- modify.movement(move_mod_3, i)
}

### 4. One leading edge ("line in the sand")
move_mod_4 <- movement
pos <- 9
move_mod_4 <- modify.movement(move_mod_4, pos)

### 5. Two leading edge ("two lines in the sand")
move_mod_5 <- movement
pos <- 9:10
for (i in pos) {
  move_mod_5 <- modify.movement(move_mod_5, i)
}

### 6. One lower, one leading edge ("squeeze / sandwich") 
move_mod_6 <- movement
pos <- c(2, 9)
for (i in pos) {
  move_mod_6 <- modify.movement(move_mod_6, i)
}

### 7. Two lower, one leading edge ("mcdouble")
move_mod_7 <- movement
pos <- c(2:3, 9)
for (i in pos) {
  move_mod_7 <- modify.movement(move_mod_7, i)
}

### 8. Two lower, two leading edge ("big mac")
move_mod_8 <- movement
pos <- c(2:3, 9:10)
for (i in pos) {
  move_mod_8 <- modify.movement(move_mod_8, i)
}