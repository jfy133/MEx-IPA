## This 'static' file of the shiny aoo is to assist in development of the shiny app
## every change made here is then transferred into the shiny app, but changing 
## variables to reactive etc. as required.

## INPUT FILE HERE



library(tidyverse)
library(gridExtra)

## input$datafile ->
## infile$datapath >- datfile
## filedata == includes ziplist, zip_filter_dirs, file_prefix, filter_list, default_runsum, last_sample, first_taxon, list_taxon_, data_damage, data_edit, data_lngt

## MAKE SELECTION HERE #############################################################
datafile <- "/home/fellows/Data/MALT-Extract_data/RII_all_off.zip"
datapath <- "/home/fellows/Data/MALT-Extract_data/RII_all_off.zip"
###################################################################################

## EDIT ##
## filedata <- reactive({
infile <- datafile ## infile <- input$datafile

  if (is.null(infile)) {
    # User has not uploaded a file yet
    return(NULL)
  }

  ## EDIT ##
  zip_list <- unzip(datapath) ## zip_list <- unzip(infile$datapath)
  zip_filter_dirs <- grep("RunSummary.txt", zip_list, value=T)
  
  file_prefix <- paste(str_split(zip_list[1], "/")[[1]][1:2], collapse="/")
  
  filter_list <- append("all", unique(str_split_fixed(zip_filter_dirs, "/", 4)[,3]))
  
  default_runsum <- suppressWarnings(
    suppressMessages(
      read_tsv(zip_list[grepl("*/default/RunSummary.txt", zip_list)])
    )
  ) %>%
    filter(Node != "Total_Count")
  list_sample <- colnames(default_runsum)[2:ncol(default_runsum)] %>% sort
  first_taxon <- colnames(default_runsum)[3]
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
  
  colnames(data_damage) <- gsub("\\(11Substitution\\)", "", colnames(data_damage))
  colnames(data_damage)[3:11] <- c("C>T_01", "C>T_02", "C>T_03", "C>T_04", 
                                   "C>T_05", "C>T_06", "C>T_07", "C>T_08",
                                   "C>T_09")
  colnames(data_damage)[23:31] <- c("D>V_01", "D>V_02", "D>V_03", "D>V_04", 
                                   "D>V_05", "D>V_06", "D>V_07", "D>V_08",
                                   "D>V_09")
  
  data_edit <- data_extractor(zip_list, "*_editDistance.txt")
  data_lngt <- data_extractor(zip_list, "*_readLengthDist.txt")
  
  read_numbers <- data_extractor(zip_list, "*RunSummary.txt")
  read_numbers <- gather(read_numbers, sample, count, 3:ncol(read_numbers)) %>% 
    rowwise %>% 
    mutate(filename = str_split_fixed(filename, "/", n = 4)[3])
  
  ## EDIT ##
  filedata <- list(default_runsum=default_runsum,
               data_damage=data_damage,
               data_edit=data_edit,
               data_lngt=data_lngt,
               file_prefix=file_prefix,
               filter_list=filter_list,
               read_numbers = read_numbers
          )
    #)
  
##})

## EDIT
# output$choose_columns1 <- renderUI({
#   # If missing input, return to avoid error later in function
#   df <- filedata()
#   
#   if(is.null(df)){
#     list_sample <- ""
#   } else {
#     # Get the data set with the appropriate name
#     list_sample <- gather(df$default_runsum,
#                           sample,
#                           count,
#                           2:ncol(df$default_runsum)) %>%
#       select(sample) %>%
#       distinct() %>%
#       arrange(sample)
#   }
#   
#   # Create the checkboxes and select them all by default
#   selectInput("sample", "Select Sample",
#               choices  = list_sample,
#               selected =list_sample)
# })

# output$choose_columns2 <- renderUI({
#   
#   # If missing input, return to avoid error later in function
#   df <- filedata()
#   
#   if(is.null(df)){
#     list_taxon <- ""
#   } else {
#     
#     # Get the data set with the appropriate name
#     list_taxon <- gather(df$default_runsum,
#                          sample,
#                          count,
#                          2:ncol(df$default_runsum)) %>%
#       select(Node) %>%
#       distinct() %>%
#       arrange(Node)
#   }
#   
#   # Create the checkboxes and select them all by default
#   
#   selectizeInput("taxon", label = "Select Taxon", choices  = list_taxon,
#                  selected = list_taxon, options = list(maxOptions = 99999))
# })
# 
# output$choose_columns3 <- renderUI({
#   
#   # If missing input, return to avoid error later in function
#   df <- filedata()
#   
#   if(is.null(df)){
#     filter_list <- ""
#   } else {
#     filter_list <- df$filter_list 
#   }
#   
#   # Create the checkboxes and select them all by default
#   
#   selectizeInput("filter", label = "Select Filter", choices  = filter_list,
#                  selected = "all", options = list(maxOptions = 99999))
# })

## MAKE SELECTIONS HERE ##########
input <- c()
input$sample <- "R32A_S0_L002_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.rma6"
input$taxon <- "Tannerella_forsythia"
input$filter <- "all"
################################
  
##output$my_plot <- #renderPlot({
  if(is.null(input$taxon)){
    NULL
  } else if(input$taxon == "") {
    NULL
  } else {
    df <- filedata#()
    
    ## Damage final data prep
    if(input$filter == "all"){
      final_damage <- df$data_damage %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Node == input$taxon)
    } else {
      final_damage <- df$data_damage %>%
        filter(grepl(input$sample, filename)) %>%
        filter(grepl(input$filter, filename)) %>%
        filter(Node == input$taxon)
  }
    

    if(input$filter == "all"){
      final_damage <- gather(final_damage,
                             Position,
                             Frequency,
                             3:ncol(final_damage)) %>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3]) %>%
        filter(dataset == "default") 
    } else {
      final_damage <- gather(final_damage,
                             Position,
                             Frequency,
                             3:ncol(final_damage)) %>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3]) %>%
        filter(dataset == input$filter)
    }
    
    final_damage <- final_damage %>%
      separate(Position, c("Modification", "Position"), "_") %>%
      filter(Modification != "considered")
    
    ### Edit Distance
    #### Sample and taxon filtering
    if(input$filter == "all"){
      final_edit <- df$data_edit %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Taxon == input$taxon)
    } else {
      final_edit <- df$data_edit %>%
        filter(grepl(input$filter, filename)) %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Taxon == input$taxon)
    }
    
    final_edit <- gather(final_edit, Distance, Reads, 3:ncol(final_edit)) %>%
      rowwise() %>%
      mutate(dataset = unlist(str_split(filename, "/"))[3])
    
    distance_levels <- c("0", "1", "2", "3", "4", "5", "6",
                         "7", "8", "9", "10", "higher")
    
    ### Read Length
    #### Sample and taxon filtering
    if(input$filter == "all"){
      final_lngt <- df$data_lngt %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Taxon == input$taxon) %>%
        gather(length, reads, 3:ncol(.)) %>%
        mutate(length = str_replace(length, "bp", "")) %>%
        mutate(length = as.numeric(length))%>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3])
    } else {
      final_lngt <- df$data_lngt %>%
        filter(grepl(input$filter, filename)) %>%
        filter(grepl(input$sample, filename)) %>%
        filter(Taxon == input$taxon) %>%
        gather(length, reads, 3:ncol(.)) %>%
        mutate(length = str_replace(length, "bp", "")) %>%
        mutate(length = as.numeric(length))%>%
        rowwise() %>%
        mutate(dataset = unlist(str_split(filename, "/"))[3])
    }
    
    ## Standard plotting info for each filter
    std_colours <- c("#F8766D", "#00B6EB")
    names(std_colours) <- c("ancient", "default")
    
    ## Plot damage
    
    fitting_data <- final_damage %>% 
      select(Modification, Position, Frequency) %>% 
      spread(Modification, Frequency) %>% 
      select(Position, `C>T`) %>% 
      filter(Position <= 10) %>% 
      mutate(Position = as.numeric(Position))
    
    fit <- nls(`C>T` ~ N * exp( - rate * Position), data = fitting_data, start = list(rate = 0.000001, N = 1), control = list(maxiter = 500))
    t <- summary(fit)$coefficients["rate",3]
    df_fit <- df.residual(fit)
    pval <- pt(t, df_fit, lower.tail = F)
    
    fit_data <- tibble(fitting_data$Position, predict(fit))
    colnames(fit_data) <- c("Position", "Fitted_point") 
    fit_data <- fit_data %>% mutate(Modification = "Predicted")
    
    print(final_damage)
    
    damage_plot <- ggplot(final_damage, aes(x=Position, y=Frequency,
                                            group=Modification,
                                            colour=Modification)) +
      geom_line(size=1.2) +
      geom_line(data = fit_data, aes(Position, Fitted_point), 
                size = 0.5,
                alpha = 0.3,
                linetype = 2,
                colour = "black") +
      theme_minimal() +
      scale_color_manual(values=c("red",
                                  "light grey",
                                  "blue",
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
            legend.position = "none") +
      ggtitle("Damage plot",
              subtitle = paste("C to T (Red), G to A (Blue)"))
    
    ## Plot edit distance
    edit_plot <- ggplot() +
      geom_bar(data=filter(final_edit, dataset == "default"), aes(factor(Distance, levels = distance_levels), Reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
      geom_bar(data=filter(final_edit, dataset == "ancient"), aes(factor(Distance, levels = distance_levels), Reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
      scale_fill_manual(values=std_colours) +
      scale_colour_manual(values=std_colours) +
      xlab("edit distance") +
      ylab("reads") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 12, face="bold")) +
      ggtitle("Edit Distance")
    
    ## Plot fragment lengths
    lngt_plot <- ggplot() +
      geom_bar(data=filter(final_lngt, dataset == "default"), aes(factor(length), reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
      geom_bar(data=filter(final_lngt, dataset == "ancient"), aes(factor(length), reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
      scale_fill_manual(values=std_colours) +
      scale_colour_manual(values=std_colours) +
      xlab("length (bp)") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 12, face="bold"),
            axis.text.x = element_text(angle = 90, hjust= 1, size = 8)) +
      ggtitle("Read Distribution")
    
    ## Plot info box
    
    read_count_no <- df$read_numbers %>%
      filter(grepl(input$sample, sample)) %>%
      filter(Node == input$taxon)
    
    final_numbers <- c()
    
    for( i in 1:nrow(read_count_no)){
      final_numbers <- append(final_numbers, 
                              paste(read_count_no[i,]$filename, 
                                    " reads: ", 
                                    read_count_no[i,]$count, 
                                    sep = ""))
    }
    
    final_numbers <- paste(final_numbers, ", ", collapse = "", sep="")
    
    
    info_sample <- paste("Displayed sample:\n", input$sample, "\n")
    info_taxon <- paste("Displayed Taxon:\n", input$taxon, "\n")
    info_numbers <- paste("Read counts per filter:\n", final_numbers, "\n")
    info_pval <- paste("Damage Exponential Fitting P-Value:\n", format(pval, scientific = FALSE), "\n")
    
    info_text <- paste(info_sample, info_taxon, info_numbers, info_pval, sep = "\n")
    
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
    
 #   outputPlot

    
    # output$down <- downloadHandler(
    #   filename =  function() {
    #     paste("malt-extract_iviewer", input$sample, "-", input$taxon,".svg")
    #   },
    #   # content is a function with argument file. content writes the plot to
    #   ## the device
    #   content = function(file) {
    #     svg(file)
    #     # draw the plot
    #     grid.arrange(damage_plot, edit_plot, lngt_plot, info_plot, nrow=2)
    #     # turn the device off
    #     dev.off()
    #   }
    # )
    
#     ## Close button
#     observe({
#       if (input$close > 0) stopApp()                             # stop shiny
#     })
#     
#   }
# ##})

  }
