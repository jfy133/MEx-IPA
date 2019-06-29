#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

######## VERSION 2 #################

library(shiny)
library(shinyFiles)
library(plotly)
library(shinycustomloader)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # Application title
    titlePanel("MEx-IPA (MALT-Extract Interactive Plotting Application)"),
    
    # Sidebar with options
    sidebarLayout(
        sidebarPanel(
            h3(strong("Options")),
            textInput("select_dir",
                      strong("Enter directory"),
                      value = NULL
            ),
            htmlOutput("report_dir"),
            uiOutput("run_button"),
            br(),
            br(),
            uiOutput("file_options"),
            uiOutput("node_options"),
            uiOutput("filter_options"),
            br(),
            textInput("remove_string",
                      strong("Remove from filename"),
                      value = NULL
            ),
            checkboxInput("interactive", 
                          "Interactive Plots?", 
                          value = FALSE, 
                          width = NULL),
            selectInput("characteristic", "Characteristic",
                        list(`DNA Damage` = "damage", 
                             `Read Length` = "length"),
                        selected = "damage")
            
            
        ),
        
        # Show a plots
        mainPanel(
            tabsetPanel(
                type = "tabs",
                tabPanel(
                    title = "Single Sample",
                    fluidRow(
                        verticalLayout(
                            br(),
                            h3("Summary Statistics"),
                            column(width = 12,
                                withLoader(DT::dataTableOutput("filterstats_plot", width = "75%"), type = "html", loader = "dnaspin"),  
                                align = "center"
                            ),
                            br(),
                            h3("Read Characteristics"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                        withLoader(uiOutput("damage_plot"), type = "html", loader = "dnaspin"),
                                        withLoader(uiOutput("length_plot"), type = "html", loader = "dnaspin")
                            ),
                            br(),
                            h3("Similarity to Reference"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                        withLoader(uiOutput("edit_plot"), type = "html", loader = "dnaspin"),
                                        withLoader(uiOutput("percentidentity_plot"), type = "html", loader = "dnaspin")
                            ),
                            br(),
                            h3("Reference Coverage"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                        withLoader(uiOutput("positionscovered_plot"), type = "html", loader = "dnaspin"),
                                        withLoader(uiOutput("coveragehist_plot"), type = "html", loader = "dnaspin")
                            )
                        )
                    )
                ),
                tabPanel(
                    title = "Multiple Samples",
                    h3("multiple"),
                    uiOutput("multisample_plots")
                )
            )
        )
    )
))
