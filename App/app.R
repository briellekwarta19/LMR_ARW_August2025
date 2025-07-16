
library(shiny)

ui <- fluidPage(
  tags$img(src = "carp_photo.png", height = "300px")
)

server <- function(input, output) {}

shinyApp(ui, server)
