######## MEx-IPA VERSION 0.4,1 #################

library(shiny)
library(DT)
library(shinycustomloader)
library(shinyWidgets)


# Define UI for application 
shinyUI(fluidPage(
    # Application title
    titlePanel("MEx-IPA (MaltExtract Interactive Plotting Application)"),
    
    # Sidebar with options
    sidebarLayout(
        sidebarPanel(
            h4("General Options"),
            textInput("select_dir",
                      strong("Enter directory"),
                      value = NULL
            ),
            htmlOutput("report_dir"),
            uiOutput("run_button"),
            br(),
            uiOutput("file_options"),
            uiOutput("node_options"),
            uiOutput("filter_options"),
            textInput("remove_string",
                      strong("Remove from filename"),
                      value = NULL
            ),
            br(),
            h4("Single Sample Options"),
            strong("Interactive"),
            switchInput("interactive", 
                          "", 
                          value = FALSE, 
                          width = "auto"),
            downloadButton('downloadPlot', 'Download PDF Report'),
            br(),
            br(),
            h4("Multiple Samples Options"),
            selectInput("characteristic", "Characteristic",
                        list(`DNA Damage` = "damage", 
                             `Read Length` = "length",
                             `Edit Distance` = "edit",
                             `Percent Identity` = "percentidentity",
                             `Positions Covered` = "positionscovered",
                             `Depth Coverage` = "coveragehist"),
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
                            h4("Summary Statistics"),
                            column(width = 12,
                                withLoader(uiOutput("filterstats_plot", width = "75%"), type = "html", loader = "dnaspin"),  
                                align = "center"
                            ),
                            br(),
                            h4("Read Characteristics"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                        withLoader(uiOutput("damage_plot"), type = "html", loader = "dnaspin"),
                                        withLoader(uiOutput("length_plot"), type = "html", loader = "dnaspin")
                            ),
                            br(),
                            h4("Similarity to Reference"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                        withLoader(uiOutput("edit_plot"), type = "html", loader = "dnaspin"),
                                        withLoader(uiOutput("percentidentity_plot"), type = "html", loader = "dnaspin")
                            ),
                            br(),
                            h4("Reference Coverage"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                        withLoader(uiOutput("positionscovered_plot"), type = "html", loader = "dnaspin"),
                                        withLoader(uiOutput("coveragehist_plot"), type = "html", loader = "dnaspin")
                            )
                        )
                    )
                ),
                tabPanel(
                    title = "Multiple Samples",
                    br(),
                    tags$b("Note"),
                    br(),
                    p("May take a few moments to load. Samples with no input data for selected taxon will not be displayed. Check single sample plot for confirmation."),
                    withLoader(uiOutput("multisample_plots"), type = "html", loader = "dnaspin")
                ),
                tabPanel(
                    title = "Multiple Taxa",
                    br(),
                    tags$b("Note"),
                    br(),
                    p("May take a few moments to load. Taxa with no input data for selected sample will not be displayed. Check single sample plot for confirmation."),
                    withLoader(uiOutput("multitaxa_plots"), type = "html", loader = "dnaspin")
                ),
                tabPanel(title = "Documentation",
                         includeMarkdown("assets/docs/documentation.md")
                )
            )
        )
    )
))
