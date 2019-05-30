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
library(shinycustomloader)

# Define UI for application that draws a histogram
shinyUI(fluidPage("MEx-IPA",
                   tabPanel("Single Sample",

    # Application title
    titlePanel("MEx-IPA (MALT-Extract Interactive Plotting Application)"),

    # Sidebar with options
    sidebarLayout(
        sidebarPanel(
            h3(strong("Options")),
            textInput("select_dir", 
                      strong("Enter directory"), 
                      value = NULL),
            br(),
            br(),
            uiOutput("file_options"),
            uiOutput("node_options"),
            uiOutput("filter_options"),
            br(),
            textInput("remove_string", 
                      strong("Remove from filename"), 
                      value = NULL)
        ),

        # Show a plots
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Single Sample",
                            h2("Input Directory"),
                            htmlOutput("report_dir"),
                            uiOutput("run_button"),
                            h2("Plots"),
                            fluidRow(
                                splitLayout(cellWidths = c("33%", "33%", "33%"), 
                                            withLoader(plotlyOutput("damage_plot"), type = "html", loader = "dnaspin"),
                                            withLoader(plotlyOutput("length_plot"), type = "html", loader = "dnaspin"), 
                                            withLoader(plotlyOutput("edit_plot"), type = "html", loader = "dnaspin")
                            ),
                            fluidRow(
                                splitLayout(cellWidths = c("33%", "33%", "33%"), 
                                            withLoader(plotlyOutput("percentidentity_plot"), type = "html", loader = "dnaspin"), 
                                            withLoader(plotlyOutput("positionscovered_plot"), type = "html", loader = "dnaspin"), 
                                            withLoader(plotlyOutput("coveragehist_plot"), type = "html", loader = "dnaspin"))
                            ),
                            br(),
                            h2("Statistics"),
                            
                            withLoader(plotOutput("filterstats_plot"), type = "html", loader = "dnaspin")
                            )
                        ),
                        tabPanel("Multiple Samples",
                                 h2("Plots"),
                                 withLoader(plotOutput("comparison_plots"), type = "html", loader = "dnaspin")
                        )
                                 
                  )
            )
        )
   )
)
)