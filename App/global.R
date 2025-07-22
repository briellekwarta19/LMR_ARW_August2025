library(shiny)
library(metapopbio)
library(bslib)

#### to do: 
#-below chunk must be run before the app. Fix this!
#add deterrents? 
#add way to add in file into the shinyapp
#multiple pages (1st page = data entry, second page = active environment)

#------------------------------------------------#
path <- here::here("LMRARW-ICPMM-AR-Demo-Files")
filename <- "/LMRARW-ICPMM-AR-Demo.xlsm"
(carp_dat <- spmm.readxl(path, filename))
## Assign objects from list
c(n_stages, n_patches, group_by, lh_order, n_timesteps, 
  stage_names, patch_names, n, matrices, MM, BB) %<-% carp_dat
