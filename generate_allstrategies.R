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
harv.vals <- seq(0,0.2,0.05)
deter_names <-  c('None', 'One lower', 'Two lower', 'One leading edge',
                  'Two leading edge', 'One lower + one leading edge ', 'Two lower + one leading edge',
                  'Two lower + two leading edge')

strategy_names <- expand.grid(harv = as.numeric_version(harv.vals), deter = deter_names)


#Create the vec-permutation matrix (dimension = stages*patches x stages*patches)
#this stays the same across all management strategies
P <- metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by)
#-------------------------------------------------------------------------------
#### Simulate strategies ####
# Harvest only strategies
A <- spmm.project.matrix(P, #vec-permutation matrix,
                         BB, #block diagonal matrix for demographic paramters 
                         MM, #block diagonal for movement paramters
                         group_by = group_by, #grouping projections (by patches here)
                         lh_order = lh_order) #order of events (demographic before movement here)

harv.vals <- seq(0,0.2,0.05)
harv.strategies <- array(0, c(length(harv.vals), n_patches))
for(i in 1:length(harv.vals)){
  harv.strategies[i,1:10] <- harv.vals[i]
}

projs <- list()

for(i in 1:length(harv.vals)){
  projs[[i]] <- spmm.project(
    n = n,  # number of stage/age animals in patch i
    A = A,
    BB = BB,
    MM = MM,
    P = P,
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches, 
    mod_mort = harv.strategies[i,]
  )
}

#### TO DO ####
#add deterrent strategies 

#### patch plot ####
#----plot for specific patch----#
#subset plot
patch_match <- rep(patch_names, each = n_stages)
ex_patch <- 'MS River & Lower L&Ds'

select <- which(strategy_names$harv == 0.1 & strategy_names$deter == 'None')

subset_harv <- projs[[select]][which(patch_match == ex_patch),]
comment(subset_harv) <- comment(projs[[select]])

spmm.plot(subset_harv, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = ex_patch)


#### all patches plot ####
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

#### calculate cost per strategy ####
calculate.cost <- function(harvest_mortality) {
  initial_cost <- 100000  # Initial cost for 0.01 harvest mortality
  growth_rate <- 1.1  # 10% increase
  cost <- ifelse(harvest_mortality == 0, 0, initial_cost * (growth_rate) ^ (harvest_mortality / 0.01))
  return(cost)
}

##### this will also be saved too ####
harv.vals <- seq(0,0.2,0.05)
deter.nums <-  c(0, 1, 2, 1, 2, 2, 3,4)

strategy_cost <- expand.grid(harv = harv.vals, deter = deter.nums)

strategy_cost <- strategy_cost %>% mutate(
  cost = calculate.cost(harv) + deter*calculate.cost(0.2)
)

##### pareto plots ####
final_pop <- final_dist <- rep(NA, length(projs))

for(i in 1:length(projs)){
  projs_select <- projs[[i]]
  final_pop[i] <- sum(projs_select[, n_timesteps])
  
  summed_pop <- projs_select[seq(1, nrow(projs_select), by = 2), ] + projs_select[seq(2, nrow(projs_select), by = 2), ]

  final_dist[i] <- sum(summed_pop[,n_timesteps] >= 0)
}

strategy_outcomes <- strategy_cost

strategy_outcomes$deter <- strategy_names$deter
strategy_outcomes$pop <- final_pop #final population
strategy_outcomes$dist <- final_dist # of patches with population

library(rPref)
test <- psel(strategy_outcomes, low(cost) * low(pop))

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
KP <- find_knee_point(strategy_outcomes$cost, strategy_outcomes$pop)

##### Save to load ######
#projs
#strategy_cost


