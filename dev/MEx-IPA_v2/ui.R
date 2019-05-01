#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking "Run App" above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFiles)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("MEx-IPA (MALT Extract Interative Plotting Application) "),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3(strong("Options")),
      br(),
      strong("Select directory"), 
      br(),
      shinyDirButton("input_dir", 
                       label = "Press here", 
                       title = "Please select a directory"),
      br(),
      br(),
      selectInput("selected_file", 
                  label = "Select sample", 
                  choices = c("sample_a", "sample_b", "sample_c"), 
                  selected = NULL, 
                  multiple = FALSE,
                  selectize = TRUE, 
                  width = NULL, 
                  size = NULL),
      selectInput("selected_node", 
                  label = "Select taxon", 
                  choices = c("taxon_a", "taxon_b", "taxon_c"), 
                  selected = NULL, 
                  multiple = FALSE,
                  selectize = TRUE, 
                  width = NULL, 
                  size = NULL),
      selectInput("selected_filter", 
                  label = "Select filter", 
                  choices  = c("ancient", "default"), 
                  selected = NULL, 
                  multiple = FALSE,
                  selectize = TRUE, 
                  width = NULL, 
                  size = NULL),
      checkboxInput("selected_interactive", 
                    label = "Interactive plots?", 
                    TRUE)
    ),
    
    # Show a plots
    mainPanel(
      h2("Plots"),
      verbatimTextOutput("input_dir", placeholder = TRUE)
    )
  )
))
