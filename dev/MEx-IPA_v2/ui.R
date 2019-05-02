#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFiles)
library(plotly)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("MEx-IPA (MALT-Extract Interactive Plotting Application"),

    # Sidebar with options
    sidebarLayout(
        sidebarPanel(
            h3(strong("Options")),
            br(),
            strong("Select directory  (not functional)"), 
            br(),
            shinyDirButton("input_dir", 
                           label = "Press here", 
                           title = "Please select a directory"),
            br(),
            br(),
            uiOutput("file_options"),
            uiOutput("node_options"),
            uiOutput("filter_options")
        ),

        # Show a plots
        mainPanel(
            
            h2("Alignment Plots"),
            fluidRow(
                splitLayout(cellWidths = c("33%", "33%", "33%"), plotlyOutput("damage_plot"), plotlyOutput("length_plot"), plotlyOutput("edit_plot")
            ),
            fluidRow(
                splitLayout(cellWidths = c("33%", "33%", "33%"), plotlyOutput("percentidentity_plot"), plotlyOutput("positionscovered_plot"), plotlyOutput("coveragehist_plot"))
            ),
            br(),
            h2("Alignment Statistics"),
            
            plotOutput("filterstats_plot")
            )
        )
    )
)
)
