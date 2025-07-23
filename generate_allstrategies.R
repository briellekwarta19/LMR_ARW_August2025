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
deter.names <-  c('None', 'One lower', 'Two lower', 'One leading edge',
                  'Two leading edge', 'One lower + one leading edge ', 'Two lower + one leading edge',
                  'Two lower + two leading edge')

strategy_names <- expand.grid(harv.vals, deter.names)
colnames(strategy_names) <- c('harv', 'deter')

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

##### Save to load ######
#projs
#strategy_cost


