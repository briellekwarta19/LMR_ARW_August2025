library("metapopbio")
library(plyr)
library(tidyverse)
library(scales)
library(gridExtra)

#Load data ---------------------------------------------------------------
path <- here::here()
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))
## Assign objects from list
c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat

filename2 <- here::here("nid-info.xlsx")
nid_info <- read.xlsx(filename2)

K <- nid_info$K

#Create data frame that lists strategy names:
harv.vals <- seq(0,0.2,0.05)
deter_names <-  c('None', 'One lower', 'Two lower', 'One leading edge',
                  'Two leading edge', 'One lower + one leading edge ', 'Two lower + one leading edge',
                  'Two lower + two leading edge')

deter_id <- seq(1:8)

strategy_names <- expand.grid(harv = harv.vals, deter = deter_id)


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
BB_harv <- list()
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
    ddf = list(K = K, beta = 0.2, theta = 1.1),
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches, 
    mod_mort = harv.strategies[i,]
  )
  
}

strategy_names[c(2,7),]


projs[[2]][,n_timesteps]
projs[[7]][,n_timesteps]

#### patch plot ####
#----plot for specific patch----#
#subset plot
patch_match <- rep(patch_names, each = n_stages)
ex_patch <- 'Nimrod Dam'

select <- which(strategy_names$harv == 0.1 & 
                  strategy_names$deter == which(deter_names == deter_names[2]))

subset_harv <- projs[[select]][which(patch_match == ex_patch),]
comment(subset_harv) <- comment(projs[[select]])

spmm.plot(subset_harv, 
          ylabs = "Rel. Abund.",
          xlabs = "Years", 
          stage_names = stage_names,
          patch_names = ex_patch)


spmm.plot(subset_harv, 
          ylabs = "Biomass",
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
final_bio <- final_dist <- rep(NA, length(projs))

for(i in 1:length(projs)){
  projs_select <- projs[[i]]
  final_bio[i] <- sum(projs_select[, n_timesteps])
  
  summed_bio <- projs_select[seq(1, nrow(projs_select), by = 2), ] + projs_select[seq(2, nrow(projs_select), by = 2), ]

  final_dist[i] <- sum(summed_bio[,n_timesteps] >= 0)
}

strategy_outcomes <- strategy_cost

strategy_outcomes$deter <- deter_names[strategy_names$deter]
strategy_outcomes$biomass <- final_bio #final population

select <- 2 #which(as.factor(strategy_outcomes$H) == input$harvs)
#select <- which(strategy_names$harv == 0.1 & 
#                  strategy_names$deter == which(deter_names == deter_names[2]))
maxbio <- 15 * 1e12 #input$maxbio * 1e12
maxcost <- 3000000 #input$maxcost* 1e6


library(rPref)
PO <- psel(strategy_outcomes, low(cost) * low(biomass))
PO 

strategy_outcomes$Strategy <- NA
strategy_outcomes$Strategy[as.numeric(rownames(PO))] <- 'Efficient Strategy'
subop <- which(strategy_outcomes$cost > maxcost | strategy_outcomes$biomass > maxbio)
strategy_outcomes$Strategy[subop] <- 'Discarded Strategy'
strategy_outcomes$Strategy[c(-(as.numeric(rownames(PO))), -subop)] <- 'Dominated Strategy'

strategy_outcomes$Strategy2 <- 'Not Selected'
strategy_outcomes$Strategy2[select] <- 'Selected Strategy'

colnames(strategy_outcomes)[5] <- "Strategy Type"
colnames(strategy_outcomes)[6] <- "Selected Strategy"

colors <- c('azure4', 'lightskyblue', 'darkgreen' )

ggplot(strategy_outcomes)+
  geom_point(aes(x = cost, y = biomass, color = `Strategy Type`, shape = `Selected Strategy`), size = 4)+
  scale_shape_manual(values=c(19, 13))+
  scale_color_manual(values = colors)+
  theme_bw() +   
  ylab("Final biomass across all patches") +
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

strategy_outcomes_table <- strategy_outcomes %>% filter(`Strategy Type` == 'Efficient Strategy' |
                                                          `Selected Strategy` == "Selected Strategy")


strategy_outcomes_table <- strategy_outcomes_table %>% arrange(`Selected Strategy`)
colnames(strategy_outcomes_table) <- c('Harvest level', 'Deterrant level', 
                                       'Cost', 'Biomass', 'Strategy Type', 'Selected Strategy')

strategy_outcomes_table$Cost <- round((strategy_outcomes_table$Cost/ 1e6),2)
strategy_outcomes_table$Biomass <- round((strategy_outcomes_table$Biomass/ 1e12),2)

colnames(strategy_outcomes_table)[3:4] <- c("Cost (M)", "Biomass (T)")

library(knitr)
library(kableExtra)


kbl <- kable(strategy_outcomes_table, "html") %>%
  kable_classic("striped", full_width = F)

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


#### Sensitivities ####
patch_nam <- rep(1:n_patches, each = n_stages)
stage_nam <- rep(c('juvenile', 'adult'), times = n_patches)

mat.names <- matrix(NA, nrow = n_patches*n_stages, ncol = n_patches*n_stages)
rownames(mat.names) <- patch_nam
colnames(mat.names) <- stage_nam
mat.names[seq(1,nrow(mat.names),by = 2),] <- 'fecund'
mat.names[seq(2,nrow(mat.names),by = 2),] <- 'surv'

mat.id <- matrix(NA, nrow = n_patches*n_stages, ncol = n_patches*n_stages)
for(i in 1:nrow(mat.id)){
  for(j in 1:ncol(mat.id)){
    mat.id[i,j] <- paste0('patch: ', rownames(mat.names)[i], ', value: ', colnames(mat.names)[j], ' ', mat.names[i,j] )
  }
}

patch_nam2r <- patch_nam2c <- rep(1:n_patches, 2)
mat.names2 <- matrix(NA, nrow = n_patches*n_stages, ncol = n_patches*n_stages)
rownames(mat.names2) <- patch_nam2r
colnames(mat.names2) <- patch_nam2c

mat.id2 <- matrix(NA, nrow = n_patches*n_stages, ncol = n_patches*n_stages)
for(r in 1:n_patches){
  for(c in 1:n_patches){
    mat.id2[r,c] <- paste0('from patch ', colnames(mat.names2)[c], ', to patch ', rownames(mat.names2)[r], ', stage: ', 'juvenile' )
  }
}

for(r in (n_patches+1):ncol(mat.id2)){
  for(c in (n_patches+1):ncol(mat.id2)){
    mat.id2[r,c] <- paste0('from patch ', colnames(mat.names2)[c], ', to patch ', rownames(mat.names2)[r], ', stage: ', 'adult' )
  }
}

BB_harv <- A_harv <- list()
testM <- testB <- rep(NA, length(strategy_names$harv))
BB_e <- top5_BB <-top5_BB_indices <- top5_BB_df <- list()
MM_e <- top5_MM <-top5_MM_indices <- top5_MM_df <- list()


for(i in 1:length(strategy_names$harv)){ 
  harv_mat <- unblk.diag(BB, n_stages)
  
  for (p in 1:10) {
    B <- harv_mat[p]
    
    M <- -log(B[[1]][-1, ])  
    M <- M + harv.strategies[i,1] #This will be adjusted per strategy
    
    B[[1]][-1, ] <- exp(-M)
    harv_mat[p] <- B
  }
  
  BB_harv[[i]] <- blk.diag(harv_mat)
  
  A_harv[[i]] <- spmm.project.matrix(P, #vec-permutation matrix,
                                BB = BB_harv[[i]], #block diagonal matrix for demographic paramters 
                                MM =  MMs[[strategy_names$deter[i]]], #block diagonal for movement paramters
                                group_by = group_by, #grouping projections (by patches here)
                                lh_order = lh_order) #order of events (demographic before movement here)
  
  #extracting elasticities for this strategy
  BB_e[[i]] <- spmm.demo.elas(BB_harv[[i]],A_harv[[i]],P, MMs[[strategy_names$deter[i]]])
  
  testB[i] <- sum(BB_e[[i]])
  
  #extracting location of top 5 rates
  top5_BB[[i]] <- order(BB_e[[i]], decreasing = TRUE)[1:5]
  #extracting indices of those rates
  top5_BB_indices[[i]] <- arrayInd(top5_BB[[i]], .dim = dim(BB_e[[i]]))
  #now combine into a dataframe which matches with the rate name:
  top5_BB_df[[i]] <- data.frame(rate = mat.id[top5_BB_indices[[i]]], value = BB_e[[i]][top5_BB[[i]]])
  top5_BB_df[[i]]$type <- i
 
  
  #do the same for MM
  MM_e[[i]] <- spmm.move.elas( MMs[[strategy_names$deter[i]]], A_harv[[i]], P, BB_harv[[i]])
  
  testM[i] <- sum(MM_e[[i]])
  
  #extracting location of top 5 rates
  top5_MM[[i]] <- order(MM_e[[i]], decreasing = TRUE)[1:5]
  #extracting indices of those rates
  top5_MM_indices[[i]] <- arrayInd(top5_MM[[i]], .dim = dim(MM_e[[i]]))
  #now combine into a dataframe which matches with the rate name:
  top5_MM_df[[i]] <- data.frame(rate = mat.id2[top5_MM_indices[[i]]], value = MM_e[[i]][top5_MM[[i]]])
  top5_MM_df[[i]]$type <- i
   
}


#combine top5 dataframes
select <- 3
top5_bb <- rbind(top5_BB_df[[1]], top5_BB_df[[select]])

# Order from base
ordered_rates <- top5_bb %>%
  filter(type == 1) %>%
  arrange(value) %>%
  pull(rate)
top5_bb$rate <- factor(top5_bb$rate, levels = ordered_rates)

top5_bb$type <- ifelse(top5_bb$type == 1, 'Base model (no management)', 'Selected strategy')

# Plot
p1 <- ggplot(top5_bb, aes(x = value, y = rate, fill = type)) +
  geom_col(position = "dodge", color = "black") +
  ggtitle("Demographic elasticities") +
  scale_fill_manual(values = c("darkslategray", "darkslategray3"))+
  labs(x = "Elasticity", y = "Demographic probability") +
  theme_minimal()+
  theme(legend.position="none")


#combine top5 dataframes
select <- 3
top5_mm <- rbind(top5_MM_df[[1]], top5_MM_df[[select]])

# Order from base
ordered_rates <- top5_mm %>%
  filter(type == 1) %>%
  arrange(value) %>%
  pull(rate)
top5_mm$rate <- factor(top5_mm$rate, levels = ordered_rates)

top5_mm$type <- ifelse(top5_mm$type == 1, 'Base model (no management)', 'Selected strategy')

# Plot
p2 <- ggplot(top5_mm, aes(x = value, y = rate, fill = type)) +
  geom_col(position = "dodge", color = "black") +
  ggtitle("Movement elasticities") +
  scale_fill_manual(values = c("darkslategray", "darkslategray3"))+
  labs(x = "Elasticity", y = "Movement probability", fill = "") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2, widths = c(0.75, 1.2))

