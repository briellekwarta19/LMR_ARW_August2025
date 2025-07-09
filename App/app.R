
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


#-----------------------------------------------#

spmm.plot2 <- function(projections, ylabs = NA, xlabs = NA, 
                      stage_names = NA, patch_names = NA) {
  comments <- comment(projections)
  group_by <- strsplit(comments, " +")[[1]][1]
  if (group_by == "patches") {
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
    
  } else if (group_by == "stages") {
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

#### UI ####
ui <- page_sidebar(
  title = "Relative Abundance Plots",
  sidebar = sidebar(
    helpText(
      "Create relative abundance plots"
    ),
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
      value = 0.1
    ),
    sliderInput(
      "deter",
      label = "Deterrant level:",
      min = 0, 
      max = 1, 
      value = 0.1
    )
  ),
  
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

#### Server ####

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Load data ---------------------------------------------------------------
  path <- here::here("LMRARW-ICPMM-AR-Demo-Files")
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
  

  #3. Make projections and plots

  output$FinalPopAll <- renderUI({
    projs <- spmm.project(
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
    HTML(paste0(
      "<b> Total final relative abundance across all patches: <b>",
      format(sum(projs[, n_timesteps]), big.mark = ",")
    ))
    
  })
  
  
  output$Patch <- renderUI({
    HTML(paste0("<b> Summary for specific patch <b>"))
  })
  
  output$Plot <- renderPlot({
    
    projs <- spmm.project(
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
  
  #add final population
  output$FinalPop <- renderUI({
    
    projs <- spmm.project(
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
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
