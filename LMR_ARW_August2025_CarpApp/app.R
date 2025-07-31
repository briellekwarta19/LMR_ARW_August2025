library(shiny)
library(bslib)
library(metapopbio)
library(dplyr)
library(tidyverse)
library(scales)
library(plyr)
library(ggrepel)
library(rPref)
library(knitr)
library(kableExtra)


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
      )[, 1])
  patch_names <-
    as.vector(
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = "metadata",
        rows = c(7:10000),
        cols = c(2),
        colNames = TRUE
      )[, 1])
  n <- as.vector(t(
    rev(openxlsx::read.xlsx(
      xlsxFile = paste0(path, filename),
      sheet = "metadata",
      rows = c(7:10000),
      cols = c(3:252),
      colNames = TRUE
    ))))
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
                      stage_names = NA, patch_names = NA,
                      ylim_max = "overall") {
  comments <- comment(projections)
  grouping <- strsplit(comments, " +")[[1]][1]
  if (grouping == "patches") {
    n_stages <- as.numeric(strsplit(comments, " +")[[2]][1])
    n_patches <- as.numeric(strsplit(comments, " +")[[3]][1])
    
    graphics::par(
      mfrow = c(min(round(n_patches / 2, 0), 3), 2),
      mar = c(5, 5, 1.5, 0.5),
      oma = rep(0.5, 4)
    )
    
    starts <- seq(1, dim(projections)[1], by = n_stages)
    ends <- c(starts - 1, dim(projections)[1])[-1]
    # throw error if starts and ends != lengths()
    for (i in seq_along(starts)) {
      idx <- starts[i]:ends[i]
      if (ylim_max == "patch") {
        ylim_vals <- c(0, round(max(projections[idx, , drop = FALSE]) + 1, -1))
      } else {
        ylim_vals <- c(0, round(max(projections) + 1, -1))
      }
      graphics::matplot(
        t(projections)[, idx],
        type = 'b',
        pch = 16,
        # col = c("black", "black"),
        ylim = ylim_vals,
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
    for (i in seq_along(starts)) {
      idx <- starts[i]:ends[i]
      if (ylim_max == "stage") {
        ylim_vals <- c(0, round(max(projections[idx, , drop = FALSE]) + 1, -1))
      } else {
        ylim_vals <- c(0, round(max(projections) + 1, -1))
      }
      graphics::matplot(
        t(projections)[, idx],
        type = 'b',
        pch = 16,
        # col = c("black", "black"),
        ylim = ylim_vals,
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
           ddf = NULL, mod_mort = NA, 
           mod_rec = NA, mod_move = NA,
           P, BB, MM) {
    
    # Initial error handling 
    lh_order <- comment(A)
    grouping <- strsplit(lh_order, " +")[[1]][1]
    A_lh_order <- strsplit(lh_order, " +")[[1]][2]
    
    try(if (is.null(comment(n)))
      stop("Please specify structure of n as either patches or stages (e.g., comment(n) <- 'patches'.')")
    )
    
    try(if (comment(n) != grouping)
      stop("Structure of n and A are not the same; both should include either 'patches' or 'stages'.")
    )
    
    # LH ORDER: patches - demo 
    if (lh_order == "patches demo") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      ## Mortality modification
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, 
                                 lh_order = A_lh_order)
      }
      ## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping,
                                 lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1])
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
      # LH ORDER: patches - move 
    } else if (lh_order == "patches move") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      ## Mortality modification
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, 
                                 lh_order = A_lh_order)
      }
      ## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping,
                                 lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1])
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
      # LH ORDER: stages - demo 
    } else if (lh_order == "stages demo") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      ## Mortality modification      
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, 
                                 lh_order = A_lh_order)
      }
      ## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping,
                                 lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1])
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
      # LH ORDER: stages - move 
    } else if (lh_order == "stages move") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      ## Mortality modificaiton      
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping, 
                                 lh_order = A_lh_order)
      }
      ## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 grouping = grouping,
                                 lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   grouping = grouping, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
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
    
    # Package output 
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
        "A projection matrix. The arg mod_move is currently ignored; please modify manually."
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

unblk.diag  <- function(blk_matrix, dimensions) {
  n_mats <- ncol(blk_matrix) / dimensions 
  
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

filename2 <- "nid-info.xlsx"
nid_info <- read.xlsx(filename2)

K <- nid_info$K

## Assign objects from list
c(n_stages, n_patches, grouping, lh_order, n_timesteps,
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

#Create data frame that lists strategy names:
harv_vals <- seq(0,0.2,0.05)

deter_names <-  c('None', 'One lower', 'Two lower', 'One leading edge',
                  'Two leading edge', 'One lower + one leading edge ', 'Two lower + one leading edge',
                  'Two lower + two leading edge')

deter_id <- seq(1:8)

strategy_names <- expand.grid(harv = harv_vals, deter = deter_id)


P <- vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          grouping = grouping)

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
                                grouping = grouping, #change to grouping
                                lh_order = lh_order)
  
  projs[[i]] <- spmm.project(
    n = n,  # number of stage/age animals in patch i
    A = A[[i]],
    BB = BB,
    MM = MMs[[strategy_names$deter[i]]],
    P = P,
    ddf = list(K = K, beta = 0.2, theta = 1.1),
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches, 
    mod_mort = harv.strategies[i,]
  )
}


# calculate cost per strategy
calculate.cost <- function(harvest_mortality) {
  initial_cost <- 100000  # Initial cost for 0.01 harvest mortality
  growth_rate <- 1.1  # 10% increase
  cost <- ifelse(harvest_mortality == 0, 0, initial_cost * (growth_rate) ^ (harvest_mortality / 0.01))
  return(cost)
}

harv_vals <- seq(0,0.2,0.05)
deter_nums <-  c(0, 1, 2, 1, 2, 2, 3,4)

strategy_cost <- expand.grid(harv = harv_vals, deter_nums = deter_nums)

strategy_cost <- strategy_cost %>% mutate(
  cost = calculate.cost(harv) + deter_nums*calculate.cost(0.2)
)

strategy_cost$cost <- ceiling(strategy_cost$cost)


# #### UI ####
ui <- navbarPage(
  title = "LMR-ARW-iCARP",
  
  ##### Introduction #####
  tabPanel("Introduction",
           fluidPage(
             h2("Invasive Carp Management App"),
             p("This app allows you to explore how different harvest and deterrent levels affect carp relative biomass across 'patches' in the Lower Mississippi River/Arkansas Red-White Rivers"),
             p("Navigate to the 'Management Strategies' tab to explore harvest and deterrent strategies"),
             p("Navigate to the 'Navigate Tradeoffs' tab to examine tradeoffs between final biomass abundance and cost across different management strategies (via Pareto efficiency analysis)"),
             p("Navigate to the 'Sensitivity Analysis' tab to explore elasticity values of key demographic and movement parameters"),
             tags$hr(),
             tags$div(
               tags$h4("Map of a subset of 'patches' in the Arkansas River", style = "text-align: center;"),
               
                # Add image here
                tags$img(src = "StudyArea.png", height = "500px", style = "display: block; margin-left: auto; margin-right: auto; margin-top: 20px; margin-bottom: 20px;"),
             ), 
             tags$hr(),
             p("Developed by Brielle Thompson & Caleb Aldridge")
           )
  ),
  
  ##### Management Strategies ####
  tabPanel("Management Strategies",
           fluidRow(
             # Sidebar (left)
             column(
               width = 3,
               wellPanel(
                 helpText("Select different harvest and deterrent levels below. Select the 'All Patches' panel to see results across all locations
                          and select the 'Single Patch' panel to see outcomes for a specific patch"),
                 sliderInput(
                   "harv",
                   label = "Harvest level:",
                   min = 0, 
                   max = 0.2, 
                   value = 0,
                   step = 0.05
                 ),
                 selectInput(
                   "deter",
                   label = "Deterrent action:",
                   choices = deter_names,
                   selected = deter_names[1]
                 )
               )
             ),
             
             # Main content (right)
             column(
               width = 9,
               fluidPage(
                 uiOutput("FinalPopAll", style = "padding-top: 20px; font-weight: bold;font-size: 23px;"),
                 tags$hr(style = "border-top: 2px solid #bbb; margin-top: 10px; margin-bottom: 20px;"),
                 uiOutput("Patch", style = "padding-top: 10px; font-weight: bold;font-size: 22px;"),
                 
                 tabsetPanel(
                   tabPanel("All Patches",
                            fluidPage(
                              tags$img(src = "VectorStudyArea.png", height = "400px", width = "100%"),
                              br(), 
                              plotOutput("AllPlot", height = "400px", width = "100%")
                            )
                   ),
                    tabPanel("Single Patch",
                             fluidPage(
                               selectInput(
                                 "patch",
                                 label = "Choose a patch to display",
                                 choices = patch_names,
                                 selected = patch_names[1]
                               ),
                               fluidRow(
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
               )
             )
           )
          )
  ),
  
  
  #### Tradeoffs #####
  #Navigate Tradeoffs
  tabPanel("Navigate Tradeoffs",
           fluidPage(
             p("Here we identify tradeoffs between two objectives: low biomass and low management cost."),
             p("To visualize this tradeoff, first select your desired strategy (the combination of harvest level and deterrant strategy) which will appear as a filled in point in the plot."),
             p("Then, you can identify your biomass constraint - referring to your maximum acceptable biomass size, and your cost constraint - referring to your maximum budget."),
             p("The resulting plot shows all management strategies. The color of the points indicate 'Strategy Type' with Discarded Strategies shown in gray (i.e., it has greater biomass and/or cost than your constraints), 
             Efficient Strategoes are in green (i.e., Pareto efficient strategies, indicating that an improvement in one objective cannot be made without worsening the other), and Dominated Strategies points are in blue (i.e., Pareto ineficient strategies,
             indicating that there is another strategy that is a least as good in one objective and better in another)"),
             p("Ideal strategies are the Efficient Strategies (i.e., the Pareto front)"),
             p("See the table below where we display the biomass and cost outcomes for the Efficient Strategies and your Selected Strategy."),
             
             tags$hr(),
             
                fluidRow(
                  column(
                    width = 3,  # Narrower sidebar
                    wellPanel(
                      helpText("Select desired harvest level and popualtion and cost constraints"),
                      
                      sliderInput(
                        "H",
                        label = "Harvest level:",
                        min = 0, 
                        max = 0.2, 
                        value = 0,
                        step = 0.05
                      ),
                      selectInput(
                        "D",
                        label = "Deterrent action:",
                        choices = deter_names,
                        selected = deter_names[1]
                      ),
                      sliderInput(
                        "maxbio",
                        label = "Biomass constraint (in trillions):",
                        min = 7, 
                        max = 15, 
                        value = 15,
                        step = 1
                      ),
                      sliderInput(
                        "maxcost",
                        label = "Cost constraint (in millions):",
                        min = 0, 
                        max = 3.5,  # Represent millions
                        value = 3,
                        step = 0.25
                      )
                    )
                  ),
                  column(
                    width = 9,
                    fluidRow(
                      column(
                        width = 9,
                        plotOutput("ParetoPlot", height = "400px", width = "100%"),
                        br(),
                        htmlOutput("ParetoTable", width = '150%')
                        )
                      )
                  )
      )
    )
  
  ),
  
  tabPanel("Sensitivity Analysis",
           fluidPage(
             h2("Add here")
           )
  )
  
  
)

#### Server ####

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ##### Final pop. and cost #####
  
  output$FinalPopAll <- renderUI({
    
    if(input$harv == 0.15){
      select <- which(strategy_names$harv == '0.15' & 
                        strategy_names$deter == which(deter_names == input$deter))
    }else{
      select <- which(strategy_names$harv == input$harv & 
                        strategy_names$deter == which(deter_names == input$deter))
    }
    

    projs_select <- projs[[select]]
    cost_select <- strategy_cost$cost[select]
    
    #Display summary:
    
    HTML(paste0(
      "<b> Total final relative biomass across all patches: <b>",
      format(sum(projs_select[, n_timesteps]), big.mark = ","), "<br>",
      "Total relative cost: $", format(cost_select, big.mark = ",", scientific = FALSE)
    ))
    
  })
  
  
  ##### Specific patch #####
  output$Patch <- renderUI({
    HTML(paste0("<b> <b>"))
  })
  
  ###### Plot ######
  output$Plot <- renderPlot({
    patch_match <- rep(patch_names, each = n_stages)
    ex_patch <- input$patch
    
    if(input$harv == 0.15){
      select <- which(strategy_names$harv == '0.15' & 
                        strategy_names$deter == which(deter_names == input$deter))
    }else{
      select <- which(strategy_names$harv == input$harv & 
                        strategy_names$deter == which(deter_names == input$deter))
    }
    
    
    subset <- projs[[select]][which(patch_match == ex_patch),]
    comment(subset) <- comment(projs[[select]])
    
    
    spmm.plot2(subset, 
               ylabs = "Rel. Abund.",
               xlabs = "Years", 
               stage_names = stage_names,
               patch_names = ex_patch)
    
  })
  
  ###### Details ######
  output$FinalPop <- renderUI({
    
    patch_match <- rep(patch_names, each = n_stages)
    ex_patch <- input$patch
    
    if(input$harv == 0.15){
      select <- which(strategy_names$harv == '0.15' & 
                        strategy_names$deter == which(deter_names == input$deter))
    }else{
      select <- which(strategy_names$harv == input$harv & 
                        strategy_names$deter == which(deter_names == input$deter))
    }
    
    
    subset <- projs[[select]][which(patch_match == ex_patch),]
    
    
    HTML(paste0(
      "<b> Final relative biomass summary</b><br>",
      "Patch selected: ", input$patch, "<br>",
      "Total biomass: ", format(sum(subset[1:2, n_timesteps]), big.mark = ","), "<br>",
      "Number of juveniles: ", format(subset[1, n_timesteps], big.mark = ","), "<br>",
      "Number of adults: ", format(subset[2, n_timesteps], big.mark = ","), "<br>"
    ))

  })
  
  
  ###### All patches plot ######
  output$AllPlot <- renderPlot({
    
    if(input$harv == 0.15){
      select <- which(strategy_names$harv == '0.15' & 
                        strategy_names$deter == which(deter_names == input$deter))
    }else{
      select <- which(strategy_names$harv == input$harv & 
                        strategy_names$deter == which(deter_names == input$deter))
    }
    
    
    projs_mat <- projs[[select]]
    
    projs_juvenile <- projs_mat[seq(1, nrow(projs_mat), by = 2), ]
    projs_adult <- projs_mat[seq(2, nrow(projs_mat), by = 2), ]
    
    patch_df_j <- adply(projs_juvenile, c(1,2))
    colnames(patch_df_j) <- c('Patch', 'Year', 'Rel.abund')
    patch_df_j$class <- 'Juvenile'
    
    patch_df_a <- adply(projs_adult, c(1,2))
    colnames(patch_df_a) <- c('Patch', 'Year', 'Rel.abund')
    patch_df_a$class <- 'Adult'
    
    patch_df <- rbind(patch_df_j, patch_df_a)
    
    ggplot(patch_df)+
      geom_point(aes(x = Year, y = Rel.abund, group = interaction(Patch, class), color = class))+
      geom_line(aes(x = Year, y = Rel.abund, group = interaction(Patch, class), color = class, linetype = class))+
      scale_color_manual(
        values = c("Juvenile" = "black", "Adult" = "salmon"),
        breaks = c("Juvenile", "Adult"),
        name = NULL)+
      scale_linetype_manual(
        values = c("Juvenile" = "solid", "Adult" = "dashed"),
        breaks = c("Juvenile", "Adult"),
        name = NULL
      ) +
      theme_bw() +   
      ylab("Relative Biomass") +
      xlab("Years")+
      theme(strip.background=element_rect(colour="white",
                                          fill="white"),
            strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
            panel.border = element_rect(colour = "gray", size = 1), 
            axis.ticks = element_blank(),
            text = element_text(size = 15)
      )+
      facet_wrap(~Patch, scales = 'free', labeller = "label_both")
    
  })
  
  ###### Pareto Plot #####
  output$ParetoPlot <- renderPlot({
    
    final_bio <- final_dist <- rep(NA, length(projs))
    
    for(i in 1:length(projs)){
      projs_select <- projs[[i]]
      final_bio[i] <- sum(projs_select[, n_timesteps])
      
      summed_bio <- projs_select[seq(1, nrow(projs_select), by = 2), ] + projs_select[seq(2, nrow(projs_select), by = 2), ]
      final_dist[i] <- sum(summed_bio[,n_timesteps] >= 0)
    }
    
    strategy_outcomes <- strategy_cost
    colnames(strategy_outcomes)[2] <- 'deter'
    
    strategy_outcomes$deter <- deter_names[strategy_names$deter]
    strategy_outcomes$biomass <- final_bio #final population
    
    
    if(input$H == 0.15){
      select <- which(strategy_names$harv == '0.15' & 
                        strategy_names$deter == which(deter_names == input$D))
    }else{
      select <- which(strategy_names$harv == input$H & 
                        strategy_names$deter == which(deter_names == input$D))
    }
    
    
    maxbio <- input$maxbio * 1e12
    maxcost <- input$maxcost* 1e6
    
    PO <- psel(strategy_outcomes, low(cost) * low(biomass))

    
    strategy_outcomes$Strategy <- NA
    strategy_outcomes$Strategy[as.numeric(rownames(PO))] <- 'Efficient Strategy'
    subop <- which(strategy_outcomes$cost > maxcost | strategy_outcomes$biomass > maxbio)
    strategy_outcomes$Strategy[subop] <- 'Discarded Strategy'
    strategy_outcomes$Strategy[c(-(as.numeric(rownames(PO))), -subop)] <- 'Dominated Strategy'
    
    strategy_outcomes$Strategy2 <- 'Not Selected'
    strategy_outcomes$Strategy2[select] <- 'Selected Strategy'
    
    colnames(strategy_outcomes)[5] <- "Strategy Type"
    colnames(strategy_outcomes)[6] <- "Selected Strategy"
    
    if(sum(strategy_outcomes$Strategy == "Discarded Strategy") > 0){
      colors <- c('azure4', 'lightskyblue', 'darkgreen' ) 
    }else{
      colors <- c('lightskyblue', 'darkgreen' )
    }
    

    
    ggplot(strategy_outcomes)+
      geom_point(aes(x = cost, y = biomass, color = `Strategy Type`, shape = `Selected Strategy`),
                 size = 4, stroke = 2)+
      scale_shape_manual(values=c(1, 19))+
      scale_color_manual(values = colors)+
      theme_bw() +   
      ylab("Final biomass across all patches") +
      ggtitle("Pareto Efficiency Plot")+
      xlab("Management cost ($)")+
      geom_hline(yintercept = maxbio, linetype = 'dashed')+
      geom_vline(xintercept = maxcost, linetype = 'dashed')+
      scale_y_continuous(labels = unit_format(unit = "T", scale = 1e-12))+
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
  
  ###### Pareto Table ######
  output$ParetoTable <- renderText({
    
    final_bio <- final_dist <- rep(NA, length(projs))
    
    for(i in 1:length(projs)){
      projs_select <- projs[[i]]
      final_bio[i] <- sum(projs_select[, n_timesteps])
      
      summed_bio <- projs_select[seq(1, nrow(projs_select), by = 2), ] + projs_select[seq(2, nrow(projs_select), by = 2), ]
      final_dist[i] <- sum(summed_bio[,n_timesteps] >= 0)
    }
    
    strategy_outcomes <- strategy_cost
    colnames(strategy_outcomes)[2] <- 'deter'
    
    strategy_outcomes$deter <- deter_names[strategy_names$deter]
    strategy_outcomes$biomass <- final_bio #final population
    
    
    if(input$H == 0.15){
      select <- which(strategy_names$harv == '0.15' & 
                        strategy_names$deter == which(deter_names == input$D))
    }else{
      select <- which(strategy_names$harv == input$H & 
                        strategy_names$deter == which(deter_names == input$D))
    }
    
    
    maxbio <- input$maxbio * 1e12
    maxcost <- input$maxcost* 1e6
    
    PO <- psel(strategy_outcomes, low(cost) * low(biomass))
    
    
    strategy_outcomes$Strategy <- NA
    strategy_outcomes$Strategy[as.numeric(rownames(PO))] <- 'Efficient Strategy'
    subop <- which(strategy_outcomes$cost > maxcost | strategy_outcomes$biomass > maxbio)
    strategy_outcomes$Strategy[subop] <- 'Discarded Strategy'
    strategy_outcomes$Strategy[c(-(as.numeric(rownames(PO))), -subop)] <- 'Dominated Strategy'
    
    strategy_outcomes$Strategy2 <- 'Not Selected'
    strategy_outcomes$Strategy2[select] <- 'Selected Strategy'
    
    colnames(strategy_outcomes)[5] <- "Strategy Type"
    colnames(strategy_outcomes)[6] <- "Selected Strategy"
    
    if(length(strategy_outcomes$Strategy == "Discarded Strategy") > 0){
      colors <- c('azure4', 'lightskyblue', 'darkgreen' ) 
    }else{
      colors <- c('lightskyblue', 'darkgreen' )
    }
    
    
    strategy_outcomes_table <- strategy_outcomes %>% filter(`Strategy Type` == 'Efficient Strategy' |
                                                              `Selected Strategy` == "Selected Strategy")
    
    
    strategy_outcomes_table <- strategy_outcomes_table %>% arrange(`Selected Strategy`)
    colnames(strategy_outcomes_table) <- c('Harvest level', 'Deterrant action', 
                                           'Cost', 'Biomass', 'Strategy Type', 'Selected Strategy')
    
    strategy_outcomes_table$Cost <- round((strategy_outcomes_table$Cost/ 1e6),2)
    strategy_outcomes_table$Biomass <- round((strategy_outcomes_table$Biomass/ 1e12),2)
    colnames(strategy_outcomes_table)[3:4] <- c("Cost (M)", "Biomass (T)")
    
    kbl <- kable(strategy_outcomes_table, "html") %>%
      #kable_classic("striped", full_width = T)
      kable_styling(full_width = F) %>% 
      column_spec(1, width_min = "12em")
    
    # Apply row coloring based on Category
    for (i in 1: length(strategy_outcomes_table$`Harvest level`)) {
      if (strategy_outcomes_table$`Strategy Type`[i] == "Discarded Strategy") {
        kbl <- kbl %>% row_spec(i, background = colors[1], bold = T)
      }
      if (strategy_outcomes_table$`Strategy Type`[i] == "Dominated Strategy") {
        kbl <- kbl %>% row_spec(i, background = colors[2], bold = T)
      }
      if (strategy_outcomes_table$`Strategy Type`[i] == "Efficient Strategy") {
        kbl <- kbl %>% row_spec(i, background = colors[3], bold = T)
      }
    }
    
    kbl
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)