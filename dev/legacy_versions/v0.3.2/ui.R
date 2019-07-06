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
                    width = NULL)
      
    ),

    # Show a plots
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Single Sample",
          fluidRow(
            verticalLayout(
                br(),
                h3("Summary Statistics"),
                withLoader(DT::dataTableOutput("filterstats_plot"), type = "html", loader = "dnaspin"),
                br(),
                h3("temp"),
                withLoader(uiOutput("damage_plot2"), type = "html", loader = "dnaspin"),
                h3("Read Characteristics"),
                splitLayout(cellWidths = c("50%", "50%"),
                            withLoader(plotlyOutput("damage_plot"), type = "html", loader = "dnaspin"),
                            withLoader(plotlyOutput("length_plot"), type = "html", loader = "dnaspin")
                            ),
                br(),
                h3("Similarity to Reference"),
                splitLayout(cellWidths = c("50%", "50%"), 
                            withLoader(plotlyOutput("edit_plot"), type = "html", loader = "dnaspin"),
                            withLoader(plotlyOutput("percentidentity_plot"), type = "html", loader = "dnaspin")
                            ),
                br(),
                h3("Reference Coverage"),
                splitLayout(cellWidths = c("50%", "50%"), 
                            withLoader(plotlyOutput("positionscovered_plot"), type = "html", loader = "dnaspin"),
                            withLoader(plotlyOutput("coveragehist_plot"), type = "html", loader = "dnaspin")
                )
            )
          )
        ),
        tabPanel(
          "Multiple Samples",
          withLoader(plotOutput("comparison_plots"), type = "html", loader = "dnaspin")
        )
      )
    )
  )
))
