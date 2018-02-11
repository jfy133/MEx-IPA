#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

options(shiny.maxRequestSize=50*1024^2)

library(shiny)
library(tidyverse)
library(gridExtra)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("MALT-Extract_iViewer"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
     sidebarPanel(
       h2("Input"),
       fileInput("datafile", "Upload ZIP File",
                 multiple = FALSE,
                 accept = c("application/zip")),
       uiOutput("choose_columns1"),
       uiOutput("choose_columns2"),
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
         "Close app")

      ),

      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("my_plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  ## Load Zipped MALT-Extract output directory and extract relevant data
  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    zip_list <- unzip(infile$datapath)

    file_prefix <- paste(str_split(zip_list[1], "/")[[1]][1:2], collapse="/")


    default_runsum <- suppressWarnings(
      suppressMessages(
        read_tsv(zip_list[grepl("*/default/RunSummary.txt", zip_list)])
        )
      ) %>%
      filter(Node != "Total_Count")
    list_sample <- colnames(default_runsum)[2:ncol(default_runsum)] %>% sort
    first_taxon <- colnames(default_runsum)[2]
    list_taxon <- select(default_runsum, Node) %>%
      distinct()

    data_extractor <- function(input, search_key){
      files_out <- zip_list[grepl(search_key, input)]
      data_out <- data_frame(filename = files_out) %>%
        mutate(file_contents = map(filename, ~ suppressWarnings(suppressMessages(read_tsv(.)))))
      data_out <- unnest(data_out)
      return(data_out)
    }

    data_damage <- data_extractor(zip_list, "*_damageMismatch")

    colnames(data_damage)[3:11] <- c("C>T_01",
                                     "C>T_02",
                                     "C>T_03",
                                     "C>T_04",
                                     "C>T_05",
                                     "C>T_06",
                                     "C>T_07",
                                     "C>T_08",
                                     "C>T_09")

    data_edit <- data_extractor(zip_list, "*_editDistance.txt")
    data_lngt <- data_extractor(zip_list, "*_readLengthDist.txt")

    return(list(default_runsum=default_runsum,
                data_damage=data_damage,
                data_edit=data_edit,
                data_lngt=data_lngt,
                file_prefix=file_prefix
                ))

  })

  output$choose_columns1 <- renderUI({
    # If missing input, return to avoid error later in function
    df <- filedata()

    if(is.null(df)){
      list_sample <- ""
    } else {
    # Get the data set with the appropriate name
    list_sample <- gather(df$default_runsum,
                          sample,
                          count,
                          2:ncol(df$default_runsum)) %>%
      select(sample) %>%
      distinct() %>%
      arrange(sample)
    }

    # Create the checkboxes and select them all by default
    selectInput("sample", "Select Sample",
                choices  = list_sample,
                selected =list_sample)
  })

  output$choose_columns2 <- renderUI({

    # If missing input, return to avoid error later in function
    df <- filedata()

    if(is.null(df)){
    list_taxon <- ""
    } else {

    # Get the data set with the appropriate name
    list_taxon <- gather(df$default_runsum,
                         sample,
                         count,
                         2:ncol(df$default_runsum)) %>%
      select(Node) %>%
      distinct() %>%
      arrange(Node)
    }

    # Create the checkboxes and select them all by default

    selectizeInput("taxon", label = "Select Taxon", choices  = list_taxon,
                   selected = list_taxon, options = list(maxOptions = 99999))
  })

  output$my_plot <- renderPlot({
    if(is.null(input$taxon)){
      NULL
    } else if(input$taxon == "") {
      NULL
    } else {
      df <- filedata()

      ## Damage final data prep
      final_damage <- df$data_damage %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Node == input$taxon)

      damage_reads <- filter(final_damage, grepl("default", filename))$considered_Matches


      final_damage <- gather(final_damage,
                             Position,
                             Frequency,
                             3:ncol(final_damage)) %>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3]) %>%
        filter(dataset == "default")

      final_damage <- as.tibble(final_damage[1:20,]) %>%
        separate(Position, c("Modification", "Position"), "_")

      ### Edit Distance
      #### Sample and taxon filtering
      final_edit <- df$data_edit %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Taxon == input$taxon)

      final_edit <- gather(final_edit, Distance, Reads, 3:ncol(final_edit)) %>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3])

      edit_reads <- final_edit %>%
        filter(dataset == "default") %>%
        summarise(sum = sum(Reads)) %>%
        sum

      distance_levels <- c("0", "1", "2", "3", "4", "5", "6",
                           "7", "8", "9", "10", "higher")

      ### Read Length
      #### Sample and taxon filtering
      final_lngt <- df$data_lngt %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Taxon == input$taxon) %>%
        gather(length, reads, 3:ncol(.)) %>%
        mutate(length = str_replace(length, "bp", "")) %>%
        mutate(length = as.numeric(length))%>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3])

      lngt_reads <- final_lngt %>%
        filter(dataset == "default") %>%
        summarise(sum = sum(reads)) %>%
        sum

      ## Plot damage
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
        xlab("position") +
        ylab("frequency") +
        theme(axis.text.x = element_text(angle = 90, hjust= 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12, face="bold"),
              legend.position = "nonw") +
        ggtitle("Damage plot",
                subtitle = paste(damage_reads,
                                 "default reads considered. C to T (Red), G to A (Blue)"))
      ## Plot edit distnace
      edit_plot <- ggplot() +
        geom_bar(data=filter(final_edit, dataset == "default"), aes(factor(Distance, levels = distance_levels), Reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
        geom_bar(data=filter(final_edit, dataset == "ancient"), aes(factor(Distance, levels = distance_levels), Reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
        xlab("edit distance") +
        ylab("reads") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12, face="bold")) +
        ggtitle("Edit Distance", subtitle = paste(edit_reads,
                                                  "default reads considered"))

      ## Plot fragment lengths
      lngt_plot <- ggplot() +
        geom_bar(data=filter(final_lngt, dataset == "default"), aes(factor(length), reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
        geom_bar(data=filter(final_lngt, dataset == "ancient"), aes(factor(length), reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
        xlab("length (bp)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12, face="bold"),
              axis.text.x = element_text(angle = 90, hjust= 1, size = 8)) +
        ggtitle("Read Distribution", subtitle = paste(lngt_reads,
                                                      "default reads considered"))
      ## Plot info box
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

      ## Combine plots and display
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

      ## Close button
      observe({
        if (input$close > 0) stopApp()                             # stop shiny
      })

    }
  })

}

# Run the application
shinyApp(ui = ui, server = server)

