library(shiny)
library(bslib)
library(metapopbio)
library(dplyr)
library(tidyverse)
library(scales)

#### to do: 
#add deterrents? 
#add way to add in file into the shinyapp

#-----------------------------------------------#
#### Functions ####
spmm.readxl <- function(path, filename) {
  sheetNames <- openxlsx::getSheetNames(paste0(path, filename))
  # metadata
  # metadata part 1
  metadata <-
    openxlsx::read.xlsx(
      xlsxFile = paste0(path, filename),
      sheet = "metadata",
      rows = c(1:5),
      cols = c(1:2),
      colNames = FALSE
    )
  for (i in 1:nrow(metadata)) {
    assign(metadata[i, 1], metadata[i, 2])
  }
  n_stages <- as.numeric(n_stages)
  n_patches <- as.numeric(n_patches)
  n_timesteps <- as.numeric(n_timesteps)
  # metadata part 2
  stage_names <-
    as.vector(
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = "metadata",
        rows = c(7:250),
        cols = c(1),
        colNames = TRUE
      )[, 1]
    )
  patch_names <-
    as.vector(
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = "metadata",
        rows = c(7:10000),
        cols = c(2),
        colNames = TRUE
      )[, 1]
    )
  n <-
    as.vector(as.matrix(
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = "metadata",
        rows = c(7:10000),
        cols = c(3:252),
        colNames = TRUE
      )[, 1]
    ))
  comment(n) <- grouping
  
  # movement matrices
  stage_idx <- grep("stage", sheetNames)
  for (i in stage_idx) {
    as.matrix(assign(
      sheetNames[i],
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = i,
        colNames = FALSE
      )
    ))
  }
  # mget(sheetNames[stage_idx])
  
  # demographic matrices
  patch_idx <- grep("patch", sheetNames)
  for (i in patch_idx) {
    assign(sheetNames[i],
           as.matrix(
             openxlsx::read.xlsx(
               xlsxFile = paste0(path, filename),
               sheet = i,
               colNames = FALSE
             )
           ))
  }
  
  B_list <- lapply(mget(sheetNames[patch_idx]), function(x) as.matrix(x))
  BB <- blk.diag(B_list)
  M_list <- lapply(mget(sheetNames[stage_idx]), function(x) as.matrix(x))
  MM <- blk.diag(M_list)
  
  return(
    list(
      n_stages = n_stages,
      n_patches = n_patches,
      grouping = grouping,
      lh_order = lh_order,
      n_timesteps = n_timesteps,
      stage_names = stage_names,
      patch_names = patch_names,
      n = n,
      mget(sheetNames[-1]),
      MM = MM,
      BB = BB
    )
  )
}

spmm.plot2 <- function(projections, ylabs = NA, xlabs = NA, 
                       stage_names = NA, patch_names = NA) {
  comments <- comment(projections)
  grouping <- strsplit(comments, " +")[[1]][1]
  if (grouping == "patches") {
    n_stages <- as.numeric(strsplit(comments, " +")[[2]][1])
    n_patches <- as.numeric(strsplit(comments, " +")[[3]][1])
    
    graphics::par(
      mfrow = c(1, 2),
      mar = c(5, 5, 1.5, 0.5),
      oma = rep(0.5, 4)
    )
    
    starts <- seq(1, dim(projections)[1], by = n_stages)
    ends <- c(starts - 1, dim(projections)[1])[-1]
    # throw error if starts and ends != lengths()
    for (i in 1:length(starts)) {
      graphics::matplot(
        t(projections)[, c(starts[i]:ends[i])],
        type = 'b',
        pch = 16,
        # col = c("black", "black"),
        ylim = c(0, round(max(projections) +1, -1)),
        ylab = ylabs,
        xlab = xlabs,
        main = paste("Patch :", patch_names[i])
      )
    }
    if (!is.na(stage_names)[1]) {
      graphics::plot.new()
      graphics::legend("center",
                       stage_names,
                       pch = 16,
                       col = 1:length(stage_names))
    }
    
  } else if (grouping == "stages") {
    n_patches <- as.numeric(strsplit(comments, " +")[[2]][1])
    n_stages <- as.numeric(strsplit(comments, " +")[[3]][1])
    
    graphics::par(
      mfrow = c(min(round(n_patches / 2, 0), 3), 2),
      mar = c(5, 5, 1.5, 0.5),
      oma = rep(0.5, 4)
    )
    
    starts <- seq(1, dim(projections)[1], by = n_patches)
    ends <- c(starts - 1, dim(projections)[1])[-1]
    # throw error if starts and ends != lengths()
    for (i in 1:length(starts)) {
      graphics::matplot(
        t(projections)[, c(starts[i]:ends[i])],
        type = 'b',
        pch = 16,
        # col = c("black", "black"),
        ylim = c(0, round(max(projections) +1 , -1)),
        ylab = ylabs,
        xlab = xlabs,
        main = paste("Stage :", stage_names[i])
      )
    }
    if (!is.na(patch_names)[1]) {
      graphics::plot.new()
      graphics::legend("center",
                       patch_names,
                       pch = 16,
                       col = 1:length(patch_names))
    }
  }
}

vec.perm <-
  function(n_stages,
           n_patches,
           grouping = c("patches", "stages")) {
    if (grouping == "patches") {
      m <- n_patches
      n <- n_stages
      P <- matrix(0, nrow = m * n, ncol = m * n)
      for (i in 1:m) {
        for (j in 1:n) {
          row_idx <- (i - 1) * n + j
          col_idx <- (j - 1) * m + i
          P[row_idx, col_idx] <- 1
        }
      }
    } else if (grouping == "stages") {
      m <- n_stages
      n <- n_patches
      P <- matrix(0, nrow = m * n, ncol = m * n)
      for (i in 1:m) {
        for (j in 1:n) {
          row_idx <- (i - 1) * n + j
          col_idx <- (j - 1) * m + i
          P[row_idx, col_idx] <- 1
        }
      }
    }
    comment(P) <- grouping
    return(P)
  }

spmm.project.matrix <- function(P, BB, MM, grouping = c("patches", "stages"), 
                                lh_order = c("demo", "move")) {
  A <- NULL  # Initialize A to NULL
  if (grouping == "patches" && lh_order == "demo") {
    A <- t(P) %*% MM %*% P %*% BB 
  } else if (grouping == "patches" && lh_order == "move") {
    A <- BB %*% t(P) %*% MM %*% P
  } else if (grouping == "stages" && lh_order == "demo") {
    A <- MM %*% P %*% BB %*% t(P)
  } else if (grouping == "stages" && lh_order == "move") {
    A <- P %*% BB %*% t(P) %*% MM
  } else {
    stop("Invalid combination of grouping and lh_order.")  # Default case
  }
  comment(A) <- paste(grouping, lh_order)
  return(A)
}

spmm.project2 <-
  function(n, A, n_timesteps,
           n_stages, n_patches, 
           ddf = NA, harv = NA, deter = NA,
           P, BB, MM) {
    
    try(if (is.null(comment(n)))
      stop("Please specify structure of n as either patches or stages (e.g., comment(n) <- 'patches'.')")
    )
    
    try(if (comment(n) != grouping)
      stop("Structure of n and A are not the same; both should include either 'patches' or 'stages'.")
    )
    
    lh_order <- comment(A)
    grouping <- strsplit(lh_order, " +")[[1]][1]
    A_lh_order <- strsplit(lh_order, " +")[[1]][2]
    
    if (lh_order == "patches demo") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (!is.na(harv)) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (length(harv) == length(M)) {
            M <- M + harv  # Add vector harv to mortality rates element-wise
          } else if (length(harv) == dim(M)[1]) {
            harv <- rep(harv, each = nrow(M))
            M <- M + harv
          } else {
            stop(
              "Length of harv must be either 1 or equal to the number of non-diagonal elements in the demographic matrix."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- dd.rec.BevertonHolt(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
            }
            if (ddf$s_type == "logistic") {
              B[[1]][1, ] <- dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "ddExponential") {
              B[[1]][1, ] <- dd.surv.exponential(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    } else if (lh_order == "patches move") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (!is.na(harv)) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (length(harv) == length(M)) {
            M <- M + harv  # Add vector harv to mortality rates element-wise
          } else if (length(harv) == dim(M)[1]) {
            harv <- rep(harv, each = nrow(M))
            M <- M + harv
          } else {
            stop(
              "Length of harv must be either 1 or equal to the number of non-diagonal elements in the demographic matrix."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              a <- B[[1]][1, ]
              B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], a, b)
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.BevertonHolt(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "logistic") {
              B[[1]][-1, ] <- B[[1]][-1, ] * dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                   grouping = grouping, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    } else if (lh_order == "stages demo") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (!is.na(harv)) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (length(harv) == length(M)) {
            M <- M + harv  # Add vector harv to mortality rates element-wise
          } else if (length(harv) == dim(M)[1]) {
            harv <- rep(harv, each = nrow(M))
            M <- M + harv
          } else {
            stop(
              "Length of harv must be either 1 or equal to the number of non-diagonal elements in the demographic matrix."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.Ricker(mat[, t - 1], ddf$r[i], ddf$K[i])
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.BevertonHolt(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "logistic") {
              B[[1]][-1, ] <- B[[1]][-1, ] * dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    } else if (lh_order == "stages move") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (!is.na(harv)) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (length(harv) == length(M)) {
            M <- M + harv  # Add vector harv to mortality rates element-wise
          } else if (length(harv) == dim(M)[1]) {
            harv <- rep(harv, each = nrow(M))
            M <- M + harv
          } else {
            stop(
              "Length of harv must be either 1 or equal to the number of non-diagonal elements in the demographic matrix."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 grouping = grouping, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.Ricker(mat[, t - 1], ddf$r[i], ddf$K[i])
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.BevertonHolt(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "logistic") {
              B[[1]][-1, ] <- B[[1]][-1, ] * dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    }
    
    if (A_lh_order == "move") {
      A_TYPE <- "movement then demography"
    } else if (A_lh_order == "demo") {
      A_TYPE <- "demography then movement"
    }
    print(
      paste(
        "Deterministic spatial matrix model projections for",
        grouping,
        "structured population vector and",
        A_TYPE,
        "A projection matrix."
      )
    )
    if (grouping == "patches") {
      comment(mat) <- paste(c(
        grouping,
        as.character(n_stages),
        as.character(n_patches)
      ))
    } else if (grouping == "stages") {
      comment(mat) <- paste(c(
        grouping,
        as.character(n_patches),
        as.character(n_stages)
      ))
    }
    return(mat)
  }

blk.diag  <- function(matlist) {
  if (!is.list(matlist)) {
    Ms <- list(...)
  } else if (is.null(matlist)) {
    stop("Please provide list of matrices!")
  }
  Ms <- matlist
  n_rows <- sum(sapply(Ms, nrow))
  n_cols <- sum(sapply(Ms, ncol))
  MM <- matrix(0, n_rows, n_cols)
  m <- 1
  n <- 1
  for (i in Ms) {
    mm <- m + nrow(i) - 1
    nn <- n + ncol(i) - 1
    MM[m:mm, n:nn] <- i
    m <- mm + 1
    n <- nn + 1
  }
  return(MM)
}

unblk.diag  <- function(blk_matrix, dimensions) {
  n_mats <- ncol(blk_matrix) / n_stages 
  
  matlist <- list()
  
  # Loop over each block and extract the submatrix
  for (i in seq_len(n_mats)) {
    start_idx <- (i - 1) * dimensions + 1
    end_idx <- start_idx + dimensions - 1
    sub_mat <- blk_matrix[start_idx:end_idx, 
                          start_idx:end_idx]
    matlist[[i]] <- sub_mat
  }
  
  return(matlist)
}

calculate.cost <- function(harvest_mortality) {
  initial_cost <- 100000  # Initial cost for 0.01 harvest mortality
  growth_rate <- 1.1  # 10% increase
  cost <- ifelse(harvest_mortality == 0, 0, initial_cost * (growth_rate) ^ (harvest_mortality / 0.01))
  return(cost)
}

find_knee_point <- function(x, y) {
  # Normalize the data
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Flip y if needed (to make both increasing)
  y_norm <- 1 - y_norm
  
  # Line from first to last point
  line_vec <- c(x_norm[length(x_norm)] - x_norm[1], y_norm[length(y_norm)] - y_norm[1])
  line_vec <- line_vec / sqrt(sum(line_vec^2))
  
  # Compute distances
  distances <- sapply(1:length(x_norm), function(i) {
    point_vec <- c(x_norm[i] - x_norm[1], y_norm[i] - y_norm[1])
    proj_len <- sum(point_vec * line_vec)
    proj_point <- line_vec * proj_len
    perp_vec <- point_vec - proj_point
    sqrt(sum(perp_vec^2))
  })
  
  
  # Find index of max distance
  knee_index <- which.max(distances)
  return(knee_index)
  
  
}

#### Data ####
filename <- "LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- metapopbio::spmm.readxl("", filename))

## Assign objects from list
c(n_stages, n_patches, grouping, lh_order, n_timesteps,
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

# #### UI ####
ui <- navbarPage(
  title = "LMR-ARW-iCARP",
  
  # Introduction Tab
  tabPanel("Introduction",
           fluidPage(
             h2("Invasive Carp Management App"),
             p("This app allows you to explore how different harvest and deterrent levels affect carp relative abundance across 'patches' in the Lower Missisippi River/Arkansas Red-White Rivers"),
             p("Navigate to the 'Explore Strategies' tab to explore harvest and deterrent strategies"),
             p("Navigate to the 'Navigate Tradeoffs' tab to identify tradeoffs between management outcomes and cost"),
             tags$hr(),
             p("Developed by Brielle Thompson & Caleb Aldridge")
           )
  ),
  
  # Explore Tab
  tabPanel("Explore Strategies",
           fluidRow(
             column(
               width = 3,  # Narrower sidebar
               wellPanel(
                 helpText("Explore outcomes of harvest and deterrent strategies"),
                 selectInput(
                   "var",
                   label = "Choose a patch to display",
                   choices = patch_names,
                   selected = patch_names[1]
                 ),
                 sliderInput(
                   "bins",
                   label = "Harvest level:",
                   min = 0, 
                   max = 1, 
                   value = 0,
                   step = 0.05
                 ),
                 sliderInput(
                   "deter",
                   label = "Deterrant level:",
                   min = 0, 
                   max = 1, 
                   value = 0, 
                   step = 0.05
                 )
               )
             ),
             column(
               width = 9,
               fluidRow(
                 column(
                   width = 12,
                   uiOutput("FinalPopAll", style = "padding-top: 20px; font-weight: bold;font-size: 23px;"),
                   tags$hr(style = "border-top: 2px solid #bbb; margin-top: 10px; margin-bottom: 20px;"),
                   uiOutput("Patch", style = "padding-top: 10px; font-weight: bold;font-size: 22px;")
                 ),
                 column(
                   width = 8,
                   plotOutput("Plot", height = "400px", width = "100%")
                 ),
                 column(
                   width = 4,
                   uiOutput("FinalPop", style = "padding-top: 20px;")
                 )
               )
             )
           )
  ), 
  # Tradeoffs Tab
  tabPanel("Navigate Tradeoffs",
           fluidPage(
             fluidRow(
               column(
                 width = 3,  # Narrower sidebar
                 wellPanel(
                   helpText("Compare chosen management strategy against other potential actions"),
                   
                   sliderInput(
                     "harvs",
                     label = "Selected harvest level:",
                     min = 0, 
                     max = 1, 
                     value = 0,
                     step = 0.05
                   ),
                   sliderInput(
                     "maxpop",
                     label = "Population constraint (maximum population):",
                     min = 0, 
                     max = 50000, 
                     value = 25000,
                     step = 1000
                   ),
                   sliderInput(
                     "maxcost",
                     label = "Cost constraint (maximum cost):",
                     min = 0, 
                     max = 1000000000, 
                     value = 500000000,
                     step = 500000
                   )
                 )
               ),
               column(
                 width = 9,
                 fluidRow(
                   column(
                     width = 8,
                     plotOutput("ParetoPlot", height = "400px", width = "100%")
                   )
                 )
               )
             )
           )
  )
  
)


#### Server ####

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ##### Final pop. and cost #####
  
  output$FinalPopAll <- renderUI({
    
    filename <- "LMRARW-ICPMM-AR-Demo.xlsm"
    (carp_dat <- metapopbio::spmm.readxl("", filename))
    c(n_stages, n_patches, grouping, lh_order, n_timesteps,
      stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat
    
    #1. creating the vec-permutation matrix (dimension = stages*patches x stages*patches)
    P <- vec.perm(n_stages = n_stages,
                  n_patches = n_patches,
                  grouping = grouping)
    
    #2. Create projection matrix
    A <- spmm.project.matrix(P, #vec-permutation matrix,
                             BB, #block diagonal matrix for demographic paramters 
                             MM, #block diagonal for movement paramters
                             grouping = grouping, #grouping projections (by patches here)
                             lh_order = lh_order) #order of events (demographic before movement here)
    
    
    
    #3. Make projections and plots
    projs <- spmm.project2(
      n = n,  # number of stage/age animals in patch i
      A = A,
      BB = BB,
      MM = MM,
      P = P,
      n_timesteps = n_timesteps,
      n_stages = n_stages,
      n_patches = n_patches,
      harv = input$bins
    )
    
    #4. Display costs:
    
    HTML(paste0(
      "<b> Total final relative abundance across all patches: <b>",
      format(sum(projs[, n_timesteps]), big.mark = ","), "<br>",
      "Total cost: $", format(ceiling(calculate.cost(input$bins)), big.mark = ",", scientific = FALSE)
    ))
    
  })
  
  
  ##### Specific patch #####
  output$Patch <- renderUI({
    HTML(paste0("<b> Summary for specific patch <b>"))
  })
  
  ###### Plot ######
  output$Plot <- renderPlot({
    
    filename <- "LMRARW-ICPMM-AR-Demo.xlsm"
    (carp_dat <- metapopbio::spmm.readxl("", filename))
    c(n_stages, n_patches, grouping, lh_order, n_timesteps,
      stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat
    
    #1. creating the vec-permutation matrix (dimension = stages*patches x stages*patches)
    P <- vec.perm(n_stages = n_stages,
                  n_patches = n_patches,
                  grouping = grouping)
    
    #2. Create projection matrix
    A <- spmm.project.matrix(P, #vec-permutation matrix,
                             BB, #block diagonal matrix for demographic paramters 
                             MM, #block diagonal for movement paramters
                             grouping = grouping, #grouping projections (by patches here)
                             lh_order = lh_order) #order of events (demographic before movement here)
    
    
    projs <- spmm.project2(
      n = n,  # number of stage/age animals in patch i
      A = A,
      BB = BB,
      MM = MM,
      P = P,
      n_timesteps = n_timesteps,
      n_stages = n_stages,
      n_patches = n_patches, 
      harv = input$bins
    )
    
    patch_match <- rep(patch_names, each = n_stages)
    ex_patch <- input$var
    
    subset <- projs[which(patch_match == ex_patch),]
    comment(subset) <- comment(projs)
    
    spmm.plot2(subset, 
               ylabs = "Rel. Abund.",
               xlabs = "Years", 
               stage_names = stage_names,
               patch_names = ex_patch)
    
  })
  
  ###### Details ######
  output$FinalPop <- renderUI({
    filename <- "LMRARW-ICPMM-AR-Demo.xlsm"
    (carp_dat <- metapopbio::spmm.readxl("", filename))
    c(n_stages, n_patches, grouping, lh_order, n_timesteps,
      stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat
    
    #1. creating the vec-permutation matrix (dimension = stages*patches x stages*patches)
    P <- vec.perm(n_stages = n_stages,
                  n_patches = n_patches,
                  grouping = grouping)
    
    #2. Create projection matrix
    A <- spmm.project.matrix(P, #vec-permutation matrix,
                             BB, #block diagonal matrix for demographic paramters 
                             MM, #block diagonal for movement paramters
                             grouping = grouping, #grouping projections (by patches here)
                             lh_order = lh_order) #order of events (demographic before movement here)
    
    
    projs <- spmm.project2(
      n = n,  # number of stage/age animals in patch i
      A = A,
      BB = BB,
      MM = MM,
      P = P,
      n_timesteps = n_timesteps,
      n_stages = n_stages,
      n_patches = n_patches,
      harv = input$bins
    )
    
    patch_match <- rep(patch_names, each = n_stages)
    ex_patch <- input$var
    
    subset <- projs[which(patch_match == ex_patch),]
    comment(subset) <- comment(projs)
    
    HTML(paste0(
      "<b> Final relative abundance summary</b><br>",
      "Patch selected: ", input$var, "<br>",
      "Total population: ", format(sum(subset[1:2, n_timesteps]), big.mark = ","), "<br>",
      "Number of juveniles: ", format(subset[1, n_timesteps], big.mark = ","), "<br>",
      "Number of adults: ", format(subset[2, n_timesteps], big.mark = ","), "<br>"
    ))
    
    
    
  })
  
  ###### Pareto Plot #####
  output$ParetoPlot <- renderPlot({
    
    filename <- "LMRARW-ICPMM-AR-Demo.xlsm"
    (carp_dat <- metapopbio::spmm.readxl("", filename))
    c(n_stages, n_patches, grouping, lh_order, n_timesteps,
      stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat
    
    #1. creating the vec-permutation matrix (dimension = stages*patches x stages*patches)
    P <- vec.perm(n_stages = n_stages,
                  n_patches = n_patches,
                  grouping = grouping)
    
    #2. Create projection matrix
    A <- spmm.project.matrix(P, #vec-permutation matrix,
                             BB, #block diagonal matrix for demographic paramters 
                             MM, #block diagonal for movement paramters
                             grouping = grouping, #grouping projections (by patches here)
                             lh_order = lh_order) #order of events (demographic before movement here)
    
    
    harvest <- seq(0,1,0.05)
    dets <- c(0)
    strategies_outcomes <- expand.grid(harvest,dets)
    colnames(strategies_outcomes) <- c('H', 'D')
    strategies_outcomes$TotalN <- NA
    strategies_outcomes$Cost <- NA
    
    for(i in 1:length(strategies_outcomes$H)){
      temp <- spmm.project2(n = n,  # number of stage/age animals in patch i
                            A = A,
                            BB = BB,
                            MM = MM,
                            P = P,
                            n_timesteps = n_timesteps,
                            n_stages = n_stages,
                            n_patches = n_patches, 
                            harv = strategies_outcomes$H[i]
      ) 
      
      strategies_outcomes$TotalN[i] <- sum(temp[,n_timesteps])
      strategies_outcomes$Cost[i] <- calculate.cost(strategies_outcomes$H[i])
      
    }
    
    KP <- find_knee_point(strategies_outcomes$Cost, strategies_outcomes$TotalN)
    
    select <- which(strategies_outcomes$H == input$harvs)
    maxpop <- input$maxpop
    maxcost <- input$maxcost
    
    strategies_outcomes$Strategy <- NA
    
    
    strategies_outcomes$Strategy[KP] <- 'Optimal Strategy'
    subop <- which(strategies_outcomes$Cost > maxcost | strategies_outcomes$TotalN > maxpop)
    strategies_outcomes$Strategy[subop] <- 'Discarded Strategy'
    strategies_outcomes$Strategy[c(-KP, -(select), -subop)] <- 'Potential Strategy'
    strategies_outcomes$Strategy[select] <- 'Selected Strategy'
    
    ggplot(strategies_outcomes)+
      geom_point(aes(x = Cost, y = TotalN, fill = Strategy), shape = 21, size = 5)+
      theme_bw() +   
      ylab("Final total population across all patches") +
      xlab("Management cost ($)")+
      geom_hline(yintercept = maxpop, linetype = 'dashed')+
      geom_vline(xintercept = maxcost, linetype = 'dashed')+
      scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6))+
      theme(strip.background=element_rect(colour="white",
                                          fill="white"),
            strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
            panel.border = element_rect(colour = "gray", size = 1.5), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 15)
      )
    
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)