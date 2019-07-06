#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

######## VERSION 0.4 #################

library(shiny)
library(tidyverse)
library(data.table)
library(plotly)
library(DT)
library(patchwork)

## Functions 
## ## Data loading

load_runsummary <- function(in_dir, result_path) {
    fread(paste0(in_dir, result_path, "RunSummary.txt")) %>% 
        as_tibble
}


load_module_files <- function(in_dir, file_ext) {
    list.files(in_dir, 
               pattern = file_ext, 
               recursive = T,
               full.names = T) %>% 
        enframe("Index", "File") %>%
        mutate(data = purrr::map(File, function(x) fread(x, header = T) %>% 
                                     as_tibble),
               File = map(File, function(file_name) {str_split(file_name, "/") %>% 
                       unlist %>% 
                       tail(n = 3) %>% 
                       paste(collapse = "/") %>% 
                       str_remove_all(file_ext)})) %>% 
        separate(File, 
                 c("Mode", "Module", "File"), 
                 sep =  "/") %>% 
        unnest %>%
        rename_at(vars(matches("Taxon")), 
                  function(x) str_replace_all(x, "Taxon", "Node")) %>%
        mutate(Mode = factor(Mode, levels = c("default", "ancient")))
}

filter_module_files <- function(in_dat, filt, remove_string, File, taxon) {
    sample <- enquo(File)
    taxon <- enquo(taxon)
    
    if (remove_string == "") {
        in_dat %>% 
            filter(File %in% !! sample,
                   Node == !! taxon,
                   Mode %in% filt)
    } else {
        in_dat %>% 
            mutate(File = str_remove(File, remove_string)) %>%
            filter(File %in% !! sample,
                   Node == !! taxon,
                   Mode %in% filt)
    }
    
}

## ## Data Cleaning functions
clean_damage <- function(x, s_filter, r_string, s_file, s_node) {
    filter_module_files(x, 
                        s_filter, 
                        r_string,
                        s_file, 
                        s_node) %>%
        gather(Mismatch, Frequency, 6:(ncol(.) - 1)) %>%
        separate(Mismatch, c("Mismatch", "Position"), "_") %>%
        mutate(Strand = if_else(Mismatch %in% c("C>T", 
                                                "D>V(11Substitution)"),
                                "5 prime", 
                                "3 prime"),
               Strand = factor(Strand, levels = c("5 prime", "3 prime")),
               Position = if_else(Position %in% names(damage_xaxis), 
                                  damage_xaxis[Position], 
                                  Position),
               Position = as_factor(Position)
        )
}

clean_length <- function(x, s_filter, r_string, s_file, s_node) {
    filter_module_files(x, 
                        s_filter, 
                        r_string,
                        s_file, 
                        s_node) %>%
        gather(Length_Bin, Alignment_Count, 6:ncol(.)) %>%
        mutate(Length_Bin = str_replace_all(Length_Bin, "bp", ""),
               Length_Bin = as.numeric(Length_Bin))
}

clean_edit <- function(x, s_filter, r_string, s_file, s_node) {
    filter_module_files(x, 
                        s_filter, 
                        r_string,
                        s_file, 
                        s_node) %>%
        gather(Edit_Distance, Alignment_Count, 6:ncol(.)) %>%
        mutate(Edit_Distance = factor(Edit_Distance, 
                                      levels = c(0:10, "higher")))
}

clean_percentidentity <- function(x, s_filter, r_string, s_file, s_node){
    filter_module_files(x, 
                        s_filter, 
                        r_string,
                        s_file, 
                        s_node) %>%
        gather(Percent_Identity, Alignment_Count, 6:ncol(.)) %>%
        mutate(Percent_Identity = factor(Percent_Identity, 
                                         levels = seq(80, 100, 5)))
}

clean_positionscovered <- function(x, s_filter, r_string, s_file, s_node){
    filter_module_files(x,
                        s_filter, 
                        r_string,
                        s_file, 
                        s_node) %>%
        select(-contains("Average"), -contains("_"), -Reference) %>%
        gather(Breadth, Percentage, 6:ncol(.)) %>%
        mutate(Breadth = str_remove_all(Breadth, "percCoveredHigher") %>% 
                   as.numeric)
}

clean_coveragehist <- function(x, s_filter, r_string, s_file, s_node){
    filter_module_files(x, 
                        s_filter, 
                        r_string,
                        s_file, 
                        s_node) %>%
        gather(Fold_Coverage, Base_Pairs, 7:ncol(.)) %>%
        mutate(Fold_Coverage = as_factor(Fold_Coverage)) %>%
        filter(Fold_Coverage != "0")
}

## ## Plotting functions

plot_damage <- function(x){
    ggplot(x, aes(Position, 
                  Frequency,
                  colour = Mismatch,
                  group = Mismatch)) +
        geom_line() + 
        labs(title = x$File %>% unique,
             subtitle = x$Node %>% unique) +
        xlab("Position (bp)") +
        ylab("Frequency of Mismatch") +
        facet_grid(File ~ Strand, scales = "free_x")  + 
        scale_colour_manual(values = mismatch_colours) +
        theme_minimal()
}

plot_col <- function(dat, xaxis, yaxis, xlabel, ylabel) {
    
    ggplot(dat, aes_string(xaxis, yaxis, fill = "Mode")) +
        geom_col() +
        labs(title = dat$File %>% unique,
             subtitle = dat$Node %>% unique) +
        xlab(xlabel) +
        ylab(ylabel) +
        facet_wrap(File ~ Node, scales = "free_x")  + 
        scale_fill_manual(values = mode_colours) + 
        theme_minimal()
}


## Default Aesthetics

mode_colours <- c(default = "#1b9e77", ancient = "#7570b3") 

mismatch_colours <- c(`C>T` = "#e41a1c", 
                      `G>A` = "#377eb8",
                      `D>V(11Substitution)` = "grey", 
                      `H>B(11Substitution)` = "grey")

damage_xaxis <- c("11" = "-10", "12" = "-9", "13" = "-8", "14" = "-7",
                  "15" = "-6", "16" = "-5", "17" = "-4", "18" = "-3",
                  "19" = "-2", "20" = "-1")

filterstats_info <- c(`Reference Length` =  "ReferenceLength", 
                      `Filtering Status` = "turnedOn?", 
                      `Downsampling Status` = "downSampling?",
                      `Unfiltered Reads (n)` = "NumberOfUnfilteredReads", 
                      `Filtered Reads (n)` = "NumberOfFilteredReads",
                      `Total Alignments (n)` = "numberOfAlignments", 
                      `Unfiltered Alignments (n)` = "NumberOfUnfilteredAlignments",
                      `Total Alignments on Reference (n)` = "TotalAlignmentsOnReference",
                      `Non-Duplicates on Reference (n)` = "nonDuplicatesonReference", 
                      `Non-Stacked (n)` = "nonStacked", 
                      `Read Distribution` = "uniquePerReference",
                      `Average Coverage on Reference (X)` = "AverageCoverage",
                      `Standard Deviation Coverage on Reference (X)` = "Coverge_StandardDeviation",
                      `Mean Read Length (bp)` = "Mean", 
                      `Standard Deviation Read Length (bp)` = "StandardDev", 
                      `Geometric Mean Read Length (bp)` = "GeometricMean", 
                      `Median Read Length (bp)` = "Median")


# Define server logic required to on the fly load/clean/plot data 
shinyServer(function(input, output) {
    
    
    ######## User Input ###########

    maltExtract_data <- eventReactive(input$submit, {
        
        input_dir <- paste0(input$select_dir,"/")
        
        if (!file.exists(paste0(input_dir, "default/RunSummary.txt"))) {
            return(NULL)
        }
        
        default_runsummary <- load_runsummary(input_dir, "default/")
        damageMismatch <- load_module_files(input_dir, "_damageMismatch.txt")
        editDistance <- load_module_files(input_dir, "_editDistance.txt")
        readLengthDist <- load_module_files(input_dir, "_readLengthDist.txt")
        filterStats <- load_module_files(input_dir, "_filterTable.txt")
        percentIdentity <- load_module_files(input_dir, "_percentIdentity.txt")
        alignmentDist <- load_module_files(input_dir, "_alignmentDist.txt")
        readLengthStat <- load_module_files(input_dir, "_readLengthStat.txt")
        coverageHist <- load_module_files(input_dir, "_coverageHist.txt")
        positionsCovered <- load_module_files(input_dir, "_postionsCovered.txt")
        
        ## Load data specific information for user options
        node_names <- default_runsummary %>% pull(Node) 
        file_names <- default_runsummary %>% colnames %>% .[. != "Node"]
        filter_names <- list("all", "default", "ancient")
        
        list(default_runsummary = default_runsummary,
             damageMismatch = damageMismatch,
             editDistance = editDistance,
             readLengthDist = readLengthDist,
             filterStats = filterStats,
             percentIdentity = percentIdentity,
             alignmentDist = alignmentDist,
             readLengthStat = readLengthStat,
             coverageHist = coverageHist,
             positionsCovered = positionsCovered,
             node_names = node_names,
             file_names = file_names,
             filter_names = filter_names)
        
    })
    

    
    output$report_dir <- renderText({
        
        input_dir <- paste0(input$select_dir,"/")
        
        if (!file.exists(paste0(input_dir, "default/RunSummary.txt"))) {
            paste("Detected maltExtract not data. Check input directory.")
        } else {
            paste("Detected maltExtract data!")
        }
        
    })
    
    output$run_button <- renderUI({
        
        input_dir <- paste0(input$select_dir,"/")
        
        if (!file.exists(paste0(input_dir, "default/RunSummary.txt"))) {
            return(NULL)
        } else {
            return(actionButton("submit", "Load data"))
        }
    })
    
    output$file_options <- renderUI({
        dat <- maltExtract_data()
        
        if (input$remove_string == "") {
            selectInput("selected_file", 
                        "Choose Sample:", 
                        as.list(dat$file_names), 
                        selectize =  T) 
        } else {
            selectInput("selected_file", 
                        "Choose Sample:", 
                        as.list(dat$file_names) %>% 
                            str_remove(., input$remove_string), 
                        selectize =  T) 
        }
    })
    
    output$node_options <- renderUI({
        dat <- maltExtract_data()
        
        selectInput("selected_node", 
                    "Choose Taxon:", 
                    as.list(dat$node_names),
                    selectize = T)
    })
    
    output$filter_options <- renderUI({
        dat <- maltExtract_data()
        selectInput("selected_filter", "Choose Filter:", as.list(dat$filter_names), 
                    selectize =  T,
                    selected = "all")
    })
    
    ######## Single Sample Plots ###########

    filterstats_data <- reactive({
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
        cat("\nSingle Sample: Loading filter stats data for", 
            input$selected_node ,
            "from", 
            input$selected_file)
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        basicstats_data <- dat$filterStats %>%
            select(-Module) %>%
            left_join(dat$readLengthStat %>% select(-Module), 
                      by = c("Index", "Mode", "File", "Node")) %>%
            left_join(dat$alignmentDist %>% select(-Module, -Reference), 
                      by = c("Index", "Mode", "File", "Node")) %>% 
            left_join(dat$positionsCovered %>% select(-Module, 
                                                      -Reference, 
                                                      -contains("perc")), 
                      by = c("Index", "Mode", "File", "Node")) %>%
            filter_module_files(selected_filter, 
                                input$remove_string,
                                input$selected_file, 
                                input$selected_node) 
        
        list(basicstats_data = basicstats_data)
    })

    damage_data <- reactive({
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      
      cat("\nSingle Sample: Loading damage data for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
      if (input$selected_filter == "all") {
        selected_filter <- c("default", "ancient")
      } else {
        selected_filter <- input$selected_filter
      }
      
      damage_data <- clean_damage(dat$damageMismatch, 
                                  selected_filter, 
                                  input$remove_string,
                                  input$selected_file,
                                  input$selected_node)
      
      list(damage_data = damage_data)
    })
    
    length_data <- reactive({
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      
      cat("\nSingle Sample: Loading length data for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
      if (input$selected_filter == "all") {
        selected_filter <- c("default", "ancient")
      } else {
        selected_filter <- input$selected_filter
      }
      
      length_data <- clean_length(dat$readLengthDist, 
                                  selected_filter, 
                                  input$remove_string,
                                  input$selected_file,
                                  input$selected_node)
      
      list(length_data = length_data)
    })
    
    edit_data <- reactive({
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      
      cat("\nSingle Sample: Loading edit distance data for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
      if (input$selected_filter == "all") {
        selected_filter <- c("default", "ancient")
      } else {
        selected_filter <- input$selected_filter
      }
      
      edit_data <- clean_edit(dat$editDistance, 
                              selected_filter, 
                              input$remove_string,
                              input$selected_file,
                              input$selected_node)
      
      list(edit_data = edit_data)
    })
    
    percentidentity_data <- reactive({
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      
      cat("\nSingle Sample: Loading percent identity data for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
      if (input$selected_filter == "all") {
        selected_filter <- c("default", "ancient")
      } else {
        selected_filter <- input$selected_filter
      }
      
      percent_data <- clean_percentidentity(dat$percentIdentity, 
                                            selected_filter, 
                                            input$remove_string,
                                            input$selected_file,
                                            input$selected_node)
      
      list(percent_data = percent_data)
    })
    
    positionscovered_data <- reactive({
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      
      cat("\nSingle Sample: Loading breadth coverage data for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
      if (input$selected_filter == "all") {
        selected_filter <- c("default", "ancient")
      } else {
        selected_filter <- input$selected_filter
      }
      
      positionscov_data <- clean_positionscovered(dat$positionsCovered, 
                                                  selected_filter, 
                                                  input$remove_string,
                                                  input$selected_file,
                                                  input$selected_node)
      
      list(positionscov_data = positionscov_data)
      
    })
    
    coveragehist_data <- reactive({
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      
        cat("\nSingle Sample: Loading depth coverage data for", 
            input$selected_node ,
            "from", 
            input$selected_file)
      
      if (input$selected_filter == "all") {
        selected_filter <- c("default", "ancient")
      } else {
        selected_filter <- input$selected_filter
      }
      
      coveragehist_data <- clean_coveragehist(dat$coverageHist, 
                                              selected_filter, 
                                              input$remove_string,
                                              input$selected_file,
                                              input$selected_node)
      
      list(coveragehist_data = coveragehist_data)
    })
    
        
    ## Single sample plotting
    
    output$filterstats_plot <- renderUI({
          
            req(filterstats_data())
      
            cat("\nSingle Sample: Plotting filter statistics plot for", 
                input$selected_node ,
                "from", 
                input$selected_file)
      
      
            basicstats_data <- filterstats_data()$basicstats_data
            
            if (nrow(basicstats_data) == 0) {
                filterstats_out <- renderText({"\nNo input data for this taxon in this sample!"})
            } else {
            
            filterstats_data <- basicstats_data %>%
                mutate(Mean = round(Mean, digits = 2),
                       GeometricMean = round(GeometricMean, digits = 2),
                       StandardDev = round(StandardDev, digits = 2)) %>%
                gather(Information, Value, 5:ncol(.)) %>%
                select(-Index) %>% 
                mutate(Information = map(Information, function(x) 
                    names(filterstats_info[filterstats_info == x])) %>% unlist) %>%
                mutate(Information = factor(Information, 
                                            levels = rev(names(filterstats_info)
                                            )
                )
                ) %>%
                arrange(Information)
            
            
            filterstats_out <-  DT::renderDataTable({filterstats_data %>% 
                select(Mode, Information, Value) %>% 
                spread(Mode, Value) %>% 
                datatable(extensions = 'Buttons', options = list(
                        dom = 'Bfrtip',
                        buttons = 
                            list('copy', list(
                                extend = 'collection',
                                filename = 'no',
                                buttons = c('csv', 'excel', 'pdf'),
                                text = 'Download'
                            ))
                    
                    )
                )
                })
            }
            
           tagList({filterstats_out})
    })

    output$damage_plot <- renderUI({
        
        req(damage_data())
      cat("\nSingle Sample: Plotting depth coverage plot for", 
          input$selected_node ,
          "from", 
          input$selected_file)
        damage_data <- damage_data()$damage_data
        
        if (any(input$selected_filter == "all")) {
            dam_plot <- plot_damage(damage_data %>% 
                                        filter(Mode == "default"))
        } else {
            dam_plot <- plot_damage(damage_data)
        }
        
        if (nrow(damage_data) == 0) {
            damage_out <- renderText({"\nNo input data for this taxon in this sample!"})
        } else if (input$interactive) {
            damage_out <- renderPlotly({dam_plot})
        } else if (!input$interactive) {
            damage_out <- renderPlot({dam_plot + theme(strip.text.y = element_blank())})
        }
        
        
        tagList({damage_out})
        
        
    })
    
 
    output$length_plot <- renderUI({
        
        req(length_data())
      
        cat("\nSingle Sample: Plotting length plot for", 
            input$selected_node ,
            "from", 
            input$selected_file)
      
        length_data <- length_data()$length_data
        
        len_plot <- plot_col(length_data, 
                               "Length_Bin", 
                               "Alignment_Count",
                               "Read Length Bins (bp)",
                               "Alignments (n)")
        
        if (nrow(length_data) == 0) {
            length_out <- renderText({"\nNo input data for this taxon in this sample!"})
        } else if (input$interactive) {
            length_out <- renderPlotly({len_plot})
        } else if (!input$interactive) {
            length_out <- renderPlot({len_plot + theme(strip.text = element_blank())})
        }
        
        tagList({length_out})
    })
    
  
    
    output$edit_plot <- renderUI({
       
        req(edit_data())
         cat("\nSingle Sample: Plotting edit distance plot for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
        edit_data <- edit_data()$edit_data
        
        edit_plot <- plot_col(edit_data, 
                              "Edit_Distance",
                              "Alignment_Count",
                              "Edit Distance",
                              "Alignments (n)")
        
        if (nrow(edit_data) == 0) {
            edit_out <- renderText({"\nNo input data for this taxon in this sample!"})
        } else if (input$interactive) {
            edit_out <- renderPlotly({edit_plot})
        } else if (!input$interactive) {
            edit_out <- renderPlot({edit_plot})
        }
        
        tagList({edit_out})
    })
    

    
    output$percentidentity_plot <- renderUI({
        
        req(percentidentity_data())
      cat("\nSingle Sample: Plotting percent identity plot for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
        percent_data <- percentidentity_data()$percent_data
        
        perc_plot <- plot_col(percent_data, 
                              "Percent_Identity",
                              "Alignment_Count",
                              "Sequence Identity (%)",
                              "Alignments (n)")
        
        if (nrow(percent_data) == 0) {
            percent_out <- renderText({"\nNo input data for this taxon in this sample!"})
        } else if (input$interactive) {
            percent_out <- renderPlotly({perc_plot})
        } else if (!input$interactive) {
            percent_out <- renderPlot({perc_plot + theme(strip.text = element_blank())})
        }
        
        tagList({percent_out})
    })
    
 
    
    output$positionscovered_plot <- renderUI({
        
        req(positionscovered_data())
      cat("\nSingle Sample: Plotting breadth coverage plot for", 
          input$selected_node ,
          "from", 
          input$selected_file)
      
        positionscov_data <- positionscovered_data()$positionscov_data

        
        positionscov_plot <- plot_col(positionscov_data, 
                                      "Breadth",
                                      "Percentage",
                                      "Fold Coverage (X)",
                                      "Percentage of Reference (%)")
        
        if (nrow( positionscov_data) == 0) {
            positionscov_out <- renderText({"\nNo input data for this taxon in this sample!"})
        } else if (input$interactive) {
            positionscov_out <- renderPlotly({positionscov_plot})
        } else if (!input$interactive) {
            positionscov_out <- renderPlot({positionscov_plot + theme(strip.text = element_blank())})
        }
        
        tagList({positionscov_out})
    })
    

    
    output$coveragehist_plot <- renderUI({
        req(coveragehist_data())
        cat("\nSingle Sample: Plotting depth coverage plot for", 
            input$selected_node ,
            "from", 
            input$selected_file)
        coveragehist_data <- coveragehist_data()$coveragehist_data
        
        coveragehist_plot <- plot_col(coveragehist_data, 
                                 "Fold_Coverage",
                                 "Base_Pairs",
                                 "Fold Coverage (X)",
                                 "Base Pairs (n)")
        
        if (nrow(coveragehist_data) == 0) {
            coveragehist_out <- renderText({"\nNo input data for this taxon in this sample!"})
        } else if (input$interactive) {
            coveragehist_out <- renderPlotly({coveragehist_plot})
        } else if (!input$interactive) {
            coveragehist_out <- renderPlot({coveragehist_plot + theme(strip.text = element_blank())})
        }
        
        tagList({coveragehist_out})
    })
    
    ######## Download Single Sample Plot ###########
    output$downloadPlot <- downloadHandler(
        filename = function() {
                    if (input$remove_string == "") {
                       file_name <- input$selected_file
                    } else {
                        file_name <- input$selected_file %>% str_remove(input$remove_string)
                    }
            
                    return(paste0("MExIPA", "-", file_name, "-", input$selected_node, '.pdf'))
                  },
        content = function(file) {
            
            nodata_message <- ggplot() + 
                annotate("text",
                         x = 4, 
                         y = 2, 
                         size = 3, 
                         label =  "\nNo input data for this taxon in this sample!") + 
                theme_void()
                
               
            
            damage_data <- damage_data()$damage_data
            length_data <- length_data()$length_data
            edit_data <- edit_data()$edit_data
            percent_data <- percentidentity_data()$percent_data
            positionscov_data <- positionscovered_data()$positionscov_data
            coveragehist_data <- coveragehist_data()$coveragehist_data
            basicstats_data <- filterstats_data()$basicstats_data
            
            if (nrow(basicstats_data) != 0) {
  
                filterstats_data <-  basicstats_data %>%
                    mutate(Mean = round(Mean, digits = 2),
                           GeometricMean = round(GeometricMean, digits = 2),
                           StandardDev = round(StandardDev, digits = 2)) %>%
                    gather(Information, Value, 5:ncol(.)) %>%
                    select(-Index) %>% 
                    mutate(Information = map(Information, function(x) 
                        names(filterstats_info[filterstats_info == x])) %>% unlist) %>%
                    mutate(Information = factor(Information, 
                                                levels = rev(names(filterstats_info)
                                                             )
                                                )) %>%
                    arrange(Information) %>%
                    select(Mode, Information, Value)
            }
            
            if (input$remove_string == "") {
                file_name <- input$selected_file
            } else {
                file_name <- input$selected_file %>% str_remove(input$remove_string)
            }
            
            title_text <- (paste(file_name, "\n", input$selected_node))
            
            title_plot <- ggplot() + 
                annotate("text", x = 4, y = 2, size = 3, label = title_text) + 
                theme_void()
            
            if (nrow(basicstats_data) != 0) {
                filterstats_data <- filterstats_data %>% 
                    pull(Information) %>% 
                    enframe(name = NULL, value = "Value") %>% 
                    mutate(Mode = "Information",
                           Value = as.character(Value)) %>% 
                    bind_rows(filterstats_data %>% 
                                  mutate(Mode = as.character(Mode))) %>% 
                    mutate(Information = if_else(is.na(Information), 
                                                 Value, as.character(Information))) %>% 
                    mutate(Mode = as_factor(Mode))
                
                filterstats_plot <- ggplot(filterstats_data, aes(x = Mode, 
                                                                 y = Information, 
                                                                 label = Value)) + 
                    geom_tile(colour = "darkgrey", fill = NA) + 
                    geom_text(size = 2.3) + 
                    scale_x_discrete(position = "top") + 
                    theme_minimal() + 
                    theme(axis.ticks = element_blank(), 
                          axis.text.y = element_blank(), 
                          panel.grid = element_blank()) + 
                    ylab("")
            }
            
            if (any(input$selected_filter == "all")) {
                dam_plot <- plot_damage(damage_data %>% 
                                            filter(Mode == "default"))
            } else {
                dam_plot <- plot_damage(damage_data  + 
                                          theme(strip.text.y = element_blank()))
            }
            
            len_plot <- plot_col(length_data, 
                                 "Length_Bin", 
                                 "Alignment_Count",
                                 "Read Length Bins (bp)",
                                 "Alignments (n)") + 
              theme(strip.text = element_blank())
            
            edit_plot <- plot_col(edit_data, 
                                  "Edit_Distance",
                                  "Alignment_Count",
                                  "Edit Distance",
                                  "Alignments (n)") + 
              theme(strip.text = element_blank())
            
            
            perc_plot <- plot_col(percent_data, 
                                  "Percent_Identity",
                                  "Alignment_Count",
                                  "Sequence Identity (%)",
                                  "Alignments (n)") + 
              theme(strip.text = element_blank())
            
            positionscov_plot <- plot_col(positionscov_data, 
                                          "Breadth",
                                          "Percentage",
                                          "Fold Coverage (X)",
                                          "Percentage of Reference (%)") + 
              theme(strip.text = element_blank())
            
            
            coveragehist_plot <- plot_col(coveragehist_data, 
                                          "Fold_Coverage",
                                          "Base_Pairs",
                                          "Fold Coverage (X)",
                                          "Base Pairs (n)") + 
              theme(strip.text = element_blank())
            
            if (nrow(basicstats_data) == 0) {
                filterstats_out <- nodata_message
            } else {
                filterstats_out <- filterstats_plot
            }
            
            if (nrow(damage_data) == 0) {
                damage_out <- nodata_message
            } else {
                damage_out <- dam_plot
            }
            
            if (nrow(length_data) == 0) {
                length_out <- nodata_message
            } else {
                length_out <- len_plot
            }
            
            if (nrow(edit_data) == 0) {
                edit_out <- nodata_message
            } else {
                edit_out <- edit_plot
            }
            
            if (nrow(percent_data) == 0) {
                percent_out <- nodata_message
            } else {
                percent_out <- perc_plot

            }
            
            if (nrow( positionscov_data) == 0) {
                positionscov_out <- nodata_message
            } else {
                positionscov_out <- positionscov_plot
            }
            
            if (nrow(coveragehist_data) == 0) {
                coveragehist_out <- nodata_message
            } else {
                coveragehist_out <- coveragehist_plot
            }
            
            report <- (title_plot) + (filterstats_out + (damage_out + length_out + 
                edit_out + percent_out + 
                positionscov_out + coveragehist_out + plot_layout(ncol = 2)) + plot_layout(ncol = 1)) + plot_layout(ncol = 1, heights = c(1,22))
            
            ggsave(file, plot = report, device = cairo_pdf, width = 250, height = 325, units = "mm")
        }
    )
    
    
    ######## Multiple Sample Plot ##########

    ## 1) Make a list of all the data for each plot
    ## 2) Make placeholders for each plot's data listed in the list
    ## 3) For each entry in list, render the plot (observe executes immediately
    ## when it detects a change)
    ## Notes: as using renderUI, sends all plots to UI at once.
    
    plotInput <- reactive({
        
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        cat("\nMultiple sample plot: Loading data for all samples with", input$selected_node)
        
        if (input$selected_filter == "all" & input$characteristic != "damage") {
            selected_filter <- c("default", "ancient")
        } else if (input$selected_filter == "all" & input$characteristic == "damage") {
            selected_filter <- "default"
        } else {
            selected_filter <- input$selected_filter
        }
        
        if (input$remove_string == "") {
            n_plot <- dat$file_names
        } else {
            n_plot <- dat$file_names %>% str_remove(input$remove_string)
        }
        
        if (input$characteristic == "damage") {
            total_data <- purrr::map(n_plot, ~ clean_damage(x = dat$damageMismatch,
                                                            s_filter = selected_filter,
                                                            r_string = input$remove_string,
                                                            s_file = .x,
                                                            s_node = input$selected_node))
            
        } else if (input$characteristic == "length") {
            total_data <- purrr::map(n_plot, ~ clean_length(x = dat$readLengthDist,
                                                            s_filter = selected_filter,
                                                            r_string = input$remove_string,
                                                            s_file = .x,
                                                            s_node = input$selected_node))
        } else if (input$characteristic == "edit") {
            total_data <- purrr::map(n_plot, ~ clean_edit(x = dat$editDistance,
                                                            s_filter = selected_filter,
                                                            r_string = input$remove_string,
                                                            s_file = .x,
                                                            s_node = input$selected_node))
        } else if (input$characteristic == "percentidentity") {
            total_data <- purrr::map(n_plot, ~ clean_percentidentity(x = dat$percentIdentity,
                                                          s_filter = selected_filter,
                                                          r_string = input$remove_string,
                                                          s_file = .x,
                                                          s_node = input$selected_node))
        } else if (input$characteristic == "positionscovered") {
            total_data <- purrr::map(n_plot, ~ clean_positionscovered(x = dat$positionsCovered,
                                                          s_filter = selected_filter,
                                                          r_string = input$remove_string,
                                                          s_file = .x,
                                                          s_node = input$selected_node))
        } else if (input$characteristic == "coveragehist") {
            total_data <- purrr::map(n_plot, ~ clean_coveragehist(x = dat$coverageHist,
                                                          s_filter = selected_filter,
                                                          r_string = input$remove_string,
                                                          s_file = .x,
                                                          s_node = input$selected_node))
        }
        
        names(total_data) <- n_plot
        
        
        ## Removes entries with empty tibbles, n_plot above will be for every 
        ## file, even if no data
        total_data <- total_data[map(total_data, nrow) > 0]
        
        ## reset n_plot
        n_plot <- names(total_data)
        
        interactive <- input$interactive
        
        
        return(list("n_plot" = n_plot, "total_data" = total_data, 
                    "interactive" = interactive))
    })
    
    
    ##### Create divs
    output$multisample_plots <- renderUI({
        
        ## make all the plot(ly) objects and place in a list
            plot_output_list <- lapply(plotInput()$n_plot, function(i) {
                cat("\nMultiple taxa plot: Allocating plot slot for", i)
                plotname <- i
                column(6, plotOutput(plotname))
            })   
        
        ## create the HTML objects that correspond to the plotly objects
        do.call(tagList, plot_output_list)
    })
    
    ## Mointor for changes in plotInput object (which is data generation above)
    observe({
            lapply(plotInput()$n_plot, function(i){
              
                cat("\nMultiple sample plot: Generating plot for", i)
                output[[i]] <- renderPlot({
                    if (input$characteristic == "damage") {
                        plot_damage(plotInput()$total_data[[i]]) + 
                        theme(strip.text.y = element_blank())
                    } else if (input$characteristic == "length") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Length_Bin", 
                                 "Alignment_Count",
                                 "Read Length Bins (bp)",
                                 "Alignments (n)") + 
                        theme(strip.text = element_blank())
                    } else if (input$characteristic == "edit") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Edit_Distance",
                                 "Alignment_Count",
                                 "Edit Distance",
                                 "Alignments (n)") + 
                        theme(strip.text = element_blank())
                    } else if (input$characteristic == "percentidentity") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Percent_Identity",
                                 "Alignment_Count",
                                 "Sequence Identity (%)",
                                 "Alignments (n)") + 
                        theme(strip.text = element_blank())
                    } else if (input$characteristic == "positionscovered") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Breadth",
                                 "Percentage",
                                 "Fold Coverage (X)",
                                 "Percentage of Reference (%)") + 
                        theme(strip.text = element_blank())
                    } else if (input$characteristic == "coveragehist") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Fold_Coverage",
                                 "Base_Pairs",
                                 "Fold Coverage (X)",
                                 "Base Pairs (n)") + 
                        theme(strip.text = element_blank())
                    }
                })
            })
        
    })
    
    ######## Multiple Taxa Plot ##########
    
    ## 1) Make a list of all the data for each plot
    ## 2) Make placeholders for each plot's data listed in the list
    ## 3) For each entry in list, render the plot (observe executes immediately
    ## when it detects a change)
    ## Notes: as using renderUI, sends all plots to UI at once.
    
    plotInput2 <- reactive({
      
      dat <- maltExtract_data()
      
      req(input$selected_filter)
      cat("\nMultiple taxa plot: Loading data for all species of ", input$selected_file)
      
      if (input$selected_filter == "all" & input$characteristic != "damage") {
        selected_filter <- c("default", "ancient")
      } else if (input$selected_filter == "all" & input$characteristic == "damage") {
        selected_filter <- "default"
      } else {
        selected_filter <- input$selected_filter
      }
      
      n_plot <- dat$node_names

      
      if (input$characteristic == "damage") {
        total_data <- purrr::map(n_plot, ~ clean_damage(x = dat$damageMismatch,
                                                        s_filter = selected_filter,
                                                        r_string = input$remove_string,
                                                        s_file = input$selected_file,
                                                        s_node = .x))
        
      } else if (input$characteristic == "length") {
        total_data <- purrr::map(n_plot, ~ clean_length(x = dat$readLengthDist,
                                                        s_filter = selected_filter,
                                                        r_string = input$remove_string,
                                                        s_file = input$selected_file,
                                                        s_node = .x))
      } else if (input$characteristic == "edit") {
        total_data <- purrr::map(n_plot, ~ clean_edit(x = dat$editDistance,
                                                      s_filter = selected_filter,
                                                      r_string = input$remove_string,
                                                      s_file = input$selected_file,
                                                      s_node = .x))
      } else if (input$characteristic == "percentidentity") {
        total_data <- purrr::map(n_plot, ~ clean_percentidentity(x = dat$percentIdentity,
                                                                 s_filter = selected_filter,
                                                                 r_string = input$remove_string,
                                                                 s_file = input$selected_file,
                                                                 s_node = .x))
      } else if (input$characteristic == "positionscovered") {
        total_data <- purrr::map(n_plot, ~ clean_positionscovered(x = dat$positionsCovered,
                                                                  s_filter = selected_filter,
                                                                  r_string = input$remove_string,
                                                                  s_file = input$selected_file,
                                                                  s_node = .x))
      } else if (input$characteristic == "coveragehist") {
        total_data <- purrr::map(n_plot, ~ clean_coveragehist(x = dat$coverageHist,
                                                              s_filter = selected_filter,
                                                              r_string = input$remove_string,
                                                              s_file = input$selected_file,
                                                              s_node = .x))
      }
      
      names(total_data) <- n_plot
      
      
      ## Removes entries with empty tibbles, n_plot above will be for every 
      ## file, even if no data
      total_data <- total_data[map(total_data, nrow) > 0]
      
      ## reset n_plot
      n_plot <- names(total_data)
      
      interactive <- input$interactive
      
      
      return(list("n_plot" = n_plot, "total_data" = total_data, 
                  "interactive" = interactive))
    })
    
    
    ##### Create divs
    output$multitaxa_plots <- renderUI({
      
      ## make all the plot(ly) objects and place in a list
      plot_output_list <- lapply(plotInput2()$n_plot, function(i) {
        cat("\nMultiple taxa plot: Allocating plot slot for", i)
        plotname <- i
        column(6, plotOutput(plotname))
      })   
      
      ## create the HTML objects that correspond to the plotly objects
      do.call(tagList, plot_output_list)
    })
    
    ## Mointor for changes in plotInput object (which is data generation above)
    observe({
      lapply(plotInput2()$n_plot, function(i){
        cat("\nMultiple taxa plot: Generating plot for", i)
        output[[i]] <- renderPlot({
          if (input$characteristic == "damage") {
            plot_damage(plotInput2()$total_data[[i]])  + 
              theme(strip.text.y = element_blank())
          } else if (input$characteristic == "length") {
            plot_col(plotInput2()$total_data[[i]],
                     "Length_Bin", 
                     "Alignment_Count",
                     "Read Length Bins (bp)",
                     "Alignments (n)")  + 
              theme(strip.text = element_blank())
          } else if (input$characteristic == "edit") {
            plot_col(plotInput2()$total_data[[i]],
                     "Edit_Distance",
                     "Alignment_Count",
                     "Edit Distance",
                     "Alignments (n)") + 
              theme(strip.text = element_blank())
          } else if (input$characteristic == "percentidentity") {
            plot_col(plotInput2()$total_data[[i]],
                     "Percent_Identity",
                     "Alignment_Count",
                     "Sequence Identity (%)",
                     "Alignments (n)") + 
              theme(strip.text = element_blank())
          } else if (input$characteristic == "positionscovered") {
            plot_col(plotInput2()$total_data[[i]],
                     "Breadth",
                     "Percentage",
                     "Fold Coverage (X)",
                     "Percentage of Reference (%)") + 
              theme(strip.text = element_blank())
          } else if (input$characteristic == "coveragehist") {
            plot_col(plotInput2()$total_data[[i]],
                     "Fold_Coverage",
                     "Base_Pairs",
                     "Fold Coverage (X)",
                     "Base Pairs (n)") + 
              theme(strip.text = element_blank())
          }
        })
      })
      
    })
    
})

