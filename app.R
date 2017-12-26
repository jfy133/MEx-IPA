#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## TO DO
## Close button - completed 2017-12-26
## Download SVG plot - completed 2017-12-26
## make data import function - pending
## make description box on how to use - completed 2017-12-26
## display general sample statistics (e.g. % of species != 0 in sample) - pending
## directory selection(?) - pending

args = commandArgs(trailingOnly=TRUE)

## Load libraries
library(shiny)
library(tidyverse)
library(gridExtra)

## Load our data directory and extract sample and taxon list

############ SELECT YOUR DIRECTORY HERE ########################################
input_dir <- "../NS_170817_malt_extract"
################################################################################


default_runsum <- read_tsv(paste(input_dir,
                                 "/default/RunSummary.txt",
                                 sep = "")) %>%
  filter(Node != "Total_Count")
list_sample <- colnames(default_runsum)[2:ncol(default_runsum)]
first_taxon <- colnames(default_runsum)[2]
list_taxon <- select(default_runsum, Node) %>%
  distinct()

## Load damage data

dir_damage <- paste(input_dir, "/default/damageMismatch", sep="")
files_damage <- dir(dir_damage, pattern = "*_damageMismatch.txt") # get file names
data_damage <- data_frame(filename = files_damage) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(dir_damage, .))))
data_damage <- unnest(data_damage)
colnames(data_damage)[3:11] <- c("C>T_01",
                                 "C>T_02",
                                 "C>T_03",
                                 "C>T_04",
                                 "C>T_05",
                                 "C>T_06",
                                 "C>T_07",
                                 "C>T_08",
                                 "C>T_09")


## Load edit distance data
dir_edit <- paste(input_dir, "/default/editDistance", sep="")
files_edit <- dir(dir_edit, pattern = "*_editDistance.txt") # get file names

data_edit <- data_frame(filename = files_edit) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(dir_edit, .))))
data_edit <- unnest(data_edit)

## Load read distribution data
dir_lngt <- paste(input_dir, "/default/readDist", sep="")
files_lngt <- dir(dir_lngt, pattern = "*_readLengthDist.txt") # get file names

data_lngt <- data_frame(filename = files_lngt) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(dir_lngt, .))))
data_lngt <- unnest(data_lngt)

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("MALT-Extract iViewer"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      h2("Input"),
      selectInput("sample", label = "Sample", list_sample),
      selectInput("taxon", label = "Taxon", list_taxon),
      p("\n"),
      p("\n"),
      h2("Output"),
      downloadButton("down", label = "Download Plot as SVG"),
      p("\n"),
      tags$button(
        id = 'close',
        type = "button",
        class = "btn action-button",
        onclick = "setTimeout(function(){window.close();},500);",  # close browser
        "Close app"
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Info",
    h2("Usage"),
     p("\n"),
     tags$ol(
       tags$li("To use, firstly modify the directory path at the beginning
       of this app file before loading to your own MALT-Extract results directory
      and save."),
      tags$li("Run the app by pressing 'Run' in the top right of code box in
      Rstudio."),
      tags$li("Select a sample and taxon of interest in the dropdown menus."),
      tags$li("To search for a particular sample or taxon, click on the
      corresponding dropdown menu and start typing the name to search."),
      tags$li("To quickly switch samples and/or species, click on the
      dropdown menu to switch, then use your keyboard's arrow keys to up or
      down, and press enter to select the next taxon.")
      ),
      p("\n"),
      h2("Reading the plots"),
      h3("Damage \n"),
      p("Describes the number of C to T and G to A nucleotides at the 5'
      end of the strand and the complementary strand, respectively.
      Fragmentation of ancient molecules leads to single stranded
      overhangs that increase the risk of demamination of cytosines,
      which are 'read' by polymerases as Thymines."),
      p("Authentic ancient
      DNA should display an exponential decrease of C to T and
      complementary G to A miscoding lesions from position 1 onwards
      into the middle of the strand."),
      h3("Edit Distance \n"),
      p("Shows the number of reads with X (1->5) number of
      nucleotide differences to the aligned reference bases. A
      decreasing number of reads from 1 to 5 edit distance values
      indicates the taxon is closely related to the reference. An
      increasing number of reads from 1 to 5 edit distance values
      indicates a taxon that is not closeely related to the reference"),
      h3("Read length"),
      p("Shows the median fragment lengths of all taxa in a sample, with the
      orange line indicating the median of the selected taxon.
      Due to hydrolytic damage over time, we expect ancient DNA molecules
      to have broken into smaller and smaller fragments, around 30-200 bp.
      Any taxon that seem unusually long outside the main distribution,
      and has low damage (see above), this could indicate possible modern
      contamination."),
      p("\n")
      ),
      tabPanel("Plots", plotOutput("my_plot"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$my_plot <- renderPlot({
    ### Damage
    #### Sample and taxon filtering
    final_damage <- data_damage %>%
      filter(filename == paste(input$sample,
                               "_damageMismatch.txt", sep="")) %>%
      filter(Node == input$taxon)
    damage_reads <- final_damage$considered_Matches
    final_damage <- gather(final_damage,
                           Position,
                           Frequency,
                           3:ncol(final_damage))
    final_damage <- as.tibble(final_damage[1:20,]) %>%
      separate(Position, c("Modification", "Position"), "_")

    ### Edit Distance
    #### Sample and taxon filtering
    final_edit <- data_edit %>%
      filter(filename == paste(input$sample, "_editDistance.txt", sep="")) %>%
      filter(Node == input$taxon)
    final_edit <- gather(final_edit, Distance, Reads, 3:ncol(final_edit))
    edit_reads <- final_edit %>% summarise(sum = sum(Reads))

    ### Read Length
    #### Sample and taxon filtering
    sample_lngt <- data_lngt %>%
      filter(filename == paste(input$sample, "_readLengthDist.txt", sep="")) %>%
      filter(Median != 0)
    no_lngt_species <- nrow(sample_lngt)
    final_lgnt <- data_lngt %>%
      filter(filename == paste(input$sample, "_readLengthDist.txt", sep="")) %>%
      filter(Node == input$taxon)

    if(nrow(final_lgnt) == 0){
      text <- paste("\n   This taxon has 0 hits in this sample")
      outputPlot <- ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_minimal() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x =element_blank(),
              axis.title.y =element_blank())

      outputPlot

    } else {
      damage_plot <- ggplot(final_damage, aes(x=Position, y=Frequency,
                                              group=Modification,
                                              colour=Modification)) +
        geom_line(size=1.2) +
        theme_bw() +
        scale_color_manual(values=c("red",
                                    "blue",
                                    "light grey",
                                    "light grey")) +
        scale_x_discrete(labels=c("01", "02", "03", "04", "05", "06", "07",
                                  "08", "09", "10", "-10", "-09", "-08", "-07",
                                  "-06", "-05", "-04", "-03", "-02", "-01")) +
        theme(axis.text.x = element_text(angle = 45, hjust= 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12, face="bold"),
              legend.position = "nonw") +
        ggtitle("Damage plot",
                subtitle = paste(damage_reads,
                                 "reads considered C to T (Red), G to A (Blue)"))

      edit_plot <- ggplot(final_edit, aes(Distance, Reads)) +
        geom_bar(stat="identity", fill="light grey") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12, face="bold")) +
        ggtitle("Edit Distance", subtitle = paste(edit_reads,
                                                  "reads considered"))

      lngt_plot <- ggplot(sample_lngt, aes(x=Median)) +
        geom_histogram(fill="light grey") +
        geom_vline(aes(xintercept=final_lgnt$Median),
                   linetype="dashed", colour="dark orange", size=1.2) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12, face="bold")) +
        ggtitle("Read Distribution", subtitle = paste("No. species total:",
                                                      nrow(sample_lngt),
                                                      "\nMedian:",
                                                      final_lgnt$Median,
                                                      "SD:",
                                                      final_lgnt$StandardDev))

      info_text <- paste("Sample:",
                         input$sample,
                         "\n\n Displayed taxon:",
                         input$taxon,
                         sep=" " )

      info_plot <- ggplot() +
        annotate("text", x = 4, y = 25, size=3, label = info_text) +
        theme_minimal() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x =element_blank(),
              axis.title.y =element_blank())

      outputPlot <- grid.arrange(damage_plot, edit_plot, lngt_plot, info_plot, nrow=2)

      outputPlot

      output$down <- downloadHandler(
        filename =  function() {
          paste("malt-extract_iviewer", input$sample, "-", input$taxon,".svg")
        },
        # content is a function with argument file. content writes the plot to
        ## the device
        content = function(file) {
          svg(file)
          # draw the plot
          grid.arrange(damage_plot, edit_plot, lngt_plot, info_plot, nrow=2)
          # turn the device off
          dev.off()
        }
      )

      observe({
        if (input$close > 0) stopApp()                             # stop shiny
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
