# Set working directory ---------------------------------------------------
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load packages -----------------------------------------------------------
if (!require("pacman")) {
  install.packages(pacman)
  library("pacman")
}
p_load(metapopbio)


# Load data ---------------------------------------------------------------
path <- here::here()
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
  #harv = rep(0.2, n_patches)
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

# spmm.plot2(subset_harv, 
#           ylabs = "Rel. Abund.",
#           xlabs = "Years", 
#           stage_names = stage_names,
#           patch_names = ex_patch)

calculate.cost <- function(harvest_mortality) {
  initial_cost <- 100000  # Initial cost for 0.01 harvest mortality
  growth_rate <- 1.1  # 10% increase
  cost <- ifelse(harvest_mortality == 0, 0, initial_cost * (growth_rate) ^ (harvest_mortality / 0.01))
  return(cost)
}

#### testing pareto plot ####
#generate a ton of outcomes: 
library(dplyr)
library(tidyverse)
library(scales)
library(rPref)

harvs <- seq(0,1,0.05)
dets <- c(0)
strategies_outcomes <- expand.grid(harvs,dets)
colnames(strategies_outcomes) <- c('H', 'D')
strategies_outcomes$TotalN <- NA
strategies_outcomes$Cost <- NA

for(i in 1:length(strategies_outcomes$H)){
   temp <- spmm.project(n = n,  # number of stage/age animals in patch i
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



test <- psel(strategies_outcomes, low(Cost) * low(TotalN))
test

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
  
#Knee point: trade-off between Cost and TotalN changes most sharply
KP <- find_knee_point(strategies_outcomes$Cost, strategies_outcomes$TotalN)

paste0('Optimal point: Harvest level = ', strategies_outcomes$H[KP],
              ', Total population = ', strategies_outcomes$TotalN[KP], 
              ', Total cost = ', strategies_outcomes$Cost[KP])

#these will vary:
select <- which(as.factor(strategies_outcomes$H) == 0.7)
maxpop <- 25000
maxcost <- 500000000

strategies_outcomes$Strategy <- NA


strategies_outcomes$Strategy[KP] <- 'Optimal Strategy'
subop <- which(strategies_outcomes$Cost > maxcost | strategies_outcomes$TotalN > maxpop)
strategies_outcomes$Strategy[subop] <- 'Discarded Strategy'
strategies_outcomes$Strategy[c(-KP, -(select), -subop)] <- 'Potential Strategy'
strategies_outcomes$Strategy[select] <- 'Selected Strategy'

ggplot(strategies_outcomes)+
  geom_point(aes(x = Cost, y = TotalN, color = Strategy))+
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


#### testing all plots ####
projs_mat <- projs

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
  ylab("Relative Abundance") +
  xlab("Years")+
  theme(strip.background=element_rect(colour="white",
                                      fill="white"),
        strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
        panel.border = element_rect(colour = "gray", size = 1), 
        axis.ticks = element_blank(),
        text = element_text(size = 15)
  )+
  facet_wrap(~Patch, scales = 'free', labeller = "label_both")




