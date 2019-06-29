#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## Libraries
library(shiny)
library(tidyverse)
library(data.table)
library(plotly)
library(DT)

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
        xlab("Position (bp)") +
        ylab("Alignments (n)") +
        facet_wrap(File ~ Strand, scales = "free_x")  + 
        scale_colour_manual(values = mismatch_colours) +
        theme_minimal()
}

plot_col <- function(dat, xaxis, yaxis, xlabel, ylabel) {
    
    ggplot(dat, aes_string(xaxis, yaxis, fill = "Mode")) +
        geom_col() +
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


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
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
    
    output$report_dir <- output$report_dir_2 <- renderText({
        
        input_dir <- paste0(input$select_dir,"/")
        
        if (!file.exists(paste0(input_dir, "default/RunSummary.txt"))) {
            paste("maltExtract data is not detected. Please check your input directory.")
        } else {
            paste("maltExtract data is detected! Press 'Load data' to view.")
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
    
    output$filterstats_plot <- DT::renderDataTable({
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        filterstats_data <- dat$filterStats %>%
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
                                input$selected_node) %>%
            mutate(Mean = round(Mean, digits = 2),
                   GeometricMean = round(Mean, digits = 2),
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
        
        filterstats_data %>% 
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

    output$damage_plot <- renderUI({
        dat <- maltExtract_data()

        req(input$selected_filter)
        
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
        
        if (any(selected_filter %in% "default")) {
            dam_plot <- plot_damage(damage_data %>% 
                                           filter(Mode == "default"))
        } else {
            dam_plot <- plot_damage(damage_data)
        }
        
        if (input$interactive) {
            damage_out <- renderPlotly({dam_plot})
        } else if (!input$interactive) {
            damage_out <- renderPlot({dam_plot})
        }
        
        tagList({damage_out})
        
    })
    
    output$length_plot <- renderUI({
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
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

        len_plot <- plot_col(length_data, 
                               "Length_Bin", 
                               "Alignment_Count",
                               "Read Length Bins (bp)",
                               "Alignments (n)")
        
        if (input$interactive) {
            length_out <- renderPlotly({len_plot})
        } else if (!input$interactive) {
            length_out <- renderPlot({len_plot})
        }
        
        tagList({length_out})
    })
    
    output$edit_plot <- renderUI({
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
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
        
        edit_plot <- plot_col(edit_data, 
                              "Edit_Distance",
                              "Alignment_Count",
                              "Edit Distance",
                              "Alignments (n)")
        
        if (input$interactive) {
            edit_out <- renderPlotly({edit_plot})
        } else if (!input$interactive) {
            edit_out <- renderPlot({edit_plot})
        }
        
        tagList({edit_out})
    })
    
    output$percentidentity_plot <- renderUI({
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
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
        
        perc_plot <- plot_col(percent_data, 
                              "Percent_Identity",
                              "Alignment_Count",
                              "Sequence Identity (%)",
                              "Alignments (n)")
        
        if (input$interactive) {
            percent_out <- renderPlotly({perc_plot})
        } else if (!input$interactive) {
            percent_out <- renderPlot({perc_plot})
        }
        
        tagList({percent_out})
    })
    
    output$positionscovered_plot <- renderUI({
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
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
        
        positionscov_plot <- plot_col(positionscov_data, 
                                      "Breadth",
                                      "Percentage",
                                      "Fold Coverage (X)",
                                      "Percentage of Reference (%)")
        
        if (input$interactive) {
            positionscov_out <- renderPlotly({positionscov_plot})
        } else if (!input$interactive) {
            positionscov_out <- renderPlot({positionscov_plot})
        }
        
        tagList({positionscov_out})
    })
    
    output$coveragehist_plot <- renderUI({
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
        
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
        
        coveragehist_plot <- plot_col(coveragehist_data, 
                                 "Fold_Coverage",
                                 "Base_Pairs",
                                 "Fold Coverage (X)",
                                 "Base Pairs (n)")
        
        if (input$interactive) {
            coveragehist_out <- renderPlotly({coveragehist_plot})
        } else if (!input$interactive) {
            coveragehist_out <- renderPlot({coveragehist_plot})
        }
        
        tagList({coveragehist_out})
    })
    
    ## 1) Make a list of all the data for each plot
    ## 2) Make placeholders for each plot's data listed in the list
    ## 3) For each entry in list, render the plot (observe executes immediately
    ## when it detects a change)
    ## Notes: as using renderUI, sends all plots to UI at once.
    
    plotInput <- reactive({
        
        dat <- maltExtract_data()
        
        req(input$selected_filter)
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        n_plot <- dat$file_names
        
        ## ERROR: pmap isn't working? so total_data not loading:
        ## x, s_filter, r_string, s_file, s_node. Map doesn't work either.
        ## Try running manually?
        
        if (input$characteristic == "damage") {
            total_data <- purrr::map(n_plot, ~ clean_damage(x = dat$damageMismatch,
                                                            s_filter = "default",
                                                            r_string = input$remove_string,
                                                            s_file = .x,
                                                            s_node = input$selected_node))
            
        } else if (input$characteristic == "length") {
            total_data <- purrr::map(n_plot, ~ clean_length(x = dat$readLengthDist,
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
        
        return(list("n_plot" = n_plot, "total_data" = total_data))
    })
    
    
    ##### Create divs######
    output$multisample_plots <- renderUI({
        
        ## make all the plot(ly) objects and place in a list
        if (input$interactive) {
            plot_output_list <- lapply(plotInput()$n_plot, function(i) {
                plotname <- i
                column(6, plotlyOutput(plotname))
            })   
        } else {
            plot_output_list <- lapply(plotInput()$n_plot, function(i) {
                plotname <- i
                column(6, plotOutput(plotname))
            })   
        }
        ## create the HTML objects that correspond to the plotly objects
        do.call(tagList, plot_output_list)
    })
    
    ## Mointor for changes in plotInput object (which is data generation above)
    observe({
        if (input$interactive) {
            lapply(plotInput()$n_plot, function(i){
                output[[i]] <- renderPlotly({
                    if (input$characteristic == "damage") {
                        plot_damage(plotInput()$total_data[[i]])
                    } else if (input$characteristic == "length") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Length_Bin", 
                                 "Alignment_Count",
                                 "Read Length Bins (bp)",
                                 "Alignments (n)")
                    }
                })
            })
        } else {
            lapply(plotInput()$n_plot, function(i){
                output[[i]] <- renderPlot({
                    if (input$characteristic == "damage") {
                        plot_damage(plotInput()$total_data[[i]])
                    } else if (input$characteristic == "length") {
                        plot_col(plotInput()$total_data[[i]],
                                 "Length_Bin", 
                                 "Alignment_Count",
                                 "Read Length Bins (bp)",
                                 "Alignments (n)")
                    }
                })
            })
        }
    })
    
})
