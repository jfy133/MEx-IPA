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

## Functions 

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

## default aesthetics

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

max_plots <- 4

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
                    selectize =  T) 
    })
    
    output$damage_plot <- renderPlotly({

        req(maltExtract_data())
        
        dat <- maltExtract_data()
    
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }

        damage_data <- filter_module_files(dat$damageMismatch, 
                                           selected_filter, 
                                           input$remove_string,
                                           input$selected_file,
                                           input$selected_node) %>%
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
        
        
        if (any(selected_filter %in% "default")) {
            damage_plot <- plot_damage(damage_data %>% 
                                           filter(Mode == "default"))
        } else {
            damage_plot <- plot_damage(damage_data)
        }
        
        damage_plot
        
    })
    
    output$length_plot <- renderPlotly({
        
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        
        length_data <- filter_module_files(dat$readLengthDist, 
                                           selected_filter, 
                                           input$remove_string,
                                           input$selected_file,
                                           input$selected_node) %>%
            gather(Length_Bin, Alignment_Count, 6:ncol(.)) %>%
            mutate(Length_Bin = str_replace_all(Length_Bin, "bp", ""),
                   Length_Bin = as.numeric(Length_Bin))
        
        length_plot <- ggplot(length_data, aes(Length_Bin, 
                                               Alignment_Count, 
                                               fill = Mode)) +
            geom_col() +
            xlab("Read Length Bins (bp)") +
            ylab("Alignments (n)") +
            facet_wrap(File ~ Node, scales = "free_x")  + 
            scale_fill_manual(values = mode_colours) + 
            theme_minimal()
        
        length_plot
    })
    
    output$edit_plot <- renderPlotly({
        
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        edit_data <- filter_module_files(dat$editDistance, 
                                         selected_filter, 
                                         input$remove_string,
                                         input$selected_file, 
                                         input$selected_node) %>%
            gather(Edit_Distance, Alignment_Count, 6:ncol(.)) %>%
            mutate(Edit_Distance = factor(Edit_Distance, 
                                          levels = c(0:10, "higher")))
        
        edit_plot <- ggplot(edit_data, 
                            aes(Edit_Distance, Alignment_Count, fill = Mode)) +
            geom_col() +
            xlab("Edit Distance") +
            ylab("Alignments (n)") +
            facet_wrap(File ~ Node, scales = "free_x")  + 
            scale_fill_manual(values = mode_colours) + 
            theme_minimal()
        
        edit_plot
    })
    
    output$percentidentity_plot <- renderPlotly({
        
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
       percentidentity_data <- filter_module_files(dat$percentIdentity, 
                            selected_filter, 
                            input$remove_string,
                            input$selected_file, 
                            input$selected_node) %>%
            gather(Percent_Identity, Alignment_Count, 6:ncol(.)) %>%
            mutate(Percent_Identity = factor(Percent_Identity, 
                                             levels = seq(80, 100, 5)))
        
        percentidentity_plot <- ggplot(percentidentity_data, 
                                       aes(Percent_Identity, 
                                           Alignment_Count, 
                                           fill = Mode)) + 
            geom_col() +
            xlab("Sequence Identity (%)") +
            ylab("Alignments (n)") +
            scale_fill_manual(values = mode_colours) + 
            facet_wrap(File ~ Node, scales = "free_x")  + 
            theme_minimal()
        
        percentidentity_plot
    })
    
    output$positionscovered_plot <- renderPlotly({ 
        
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        positionscovered_data <- filter_module_files(dat$positionsCovered,
                                                     selected_filter, 
                                                     input$remove_string,
                                                     input$selected_file, 
                                                     input$selected_node) %>%
            select(-contains("Average"), -contains("_"), -Reference) %>%
            gather(Breadth, Percentage, 6:ncol(.)) %>%
            mutate(Breadth = str_remove_all(Breadth, "percCoveredHigher") %>% 
                       as.numeric)
        
        positionscovered_plot <- ggplot(positionscovered_data, 
                                        aes(Breadth, Percentage, 
                                            fill = Mode)) +
            geom_col() +
            xlab("Fold Coverage (X)") +
            ylab("Percentage of Reference (%)") +
            scale_fill_manual(values = mode_colours) + 
            facet_wrap(File ~ Node, scales = "free_x")  + 
            theme_minimal()
        
        positionscovered_plot
        
        })
        
    
    output$coveragehist_plot <- renderPlotly({ 
        
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        if (input$selected_filter == "all") {
            selected_filter <- c("default", "ancient")
        } else {
            selected_filter <- input$selected_filter
        }
        
        coveragehist_data <- filter_module_files(dat$coverageHist, 
                                                 selected_filter, 
                                                 input$remove_string,
                                                 input$selected_file, 
                                                 input$selected_node) %>%
            gather(Fold_Coverage, Base_Pairs, 7:ncol(.)) %>%
            mutate(Fold_Coverage = as_factor(Fold_Coverage)) %>%
            filter(Fold_Coverage != "0")
        
        coveragehist_plot <- ggplot(coveragehist_data, 
                                    aes(Fold_Coverage, 
                                        Base_Pairs, 
                                        fill = Mode)) +
            geom_col() +
            xlab("Fold Coverage (X)") +
            ylab("Base Pairs (n)") +
            facet_wrap(File ~ Node, scales = "free_x")  + 
            scale_fill_manual(values = mode_colours) + 
            theme_minimal()
        
        coveragehist_plot
        
        })
        
    
    output$filterstats_plot <- renderPlot({
        
        req(maltExtract_data())
        
        dat <- maltExtract_data()
        
        
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
        
        
        filterstats_plot <- ggplot(filterstats_data, 
                                   aes(y = Information,
                                       x = Mode, 
                                       label = Value)) +
            geom_tile(colour = "darkgrey", fill = NA, size = 0.3) +
            geom_text() +
            scale_x_discrete(position = "top") +
            theme_minimal() +
            theme(panel.grid = element_blank())
        
        filterstats_plot
    })
    
        output$comparison_plots <- renderPlot({
            
            req(maltExtract_data())
            
            dat <- maltExtract_data()
            
            if (input$selected_filter == "all") {
                selected_filter <- c("default", "ancient")
            } else {
                selected_filter <- input$selected_filter
            }
            
            if (input$remove_string == "") {
                file_names <- dat$file_names
            } else {
                file_names <- dat$file_names %>% 
                    str_remove(., input$remove_string)
            }

            damage_data_comparison <- filter_module_files(dat$damageMismatch, 
                                               selected_filter, 
                                               input$remove_string,
                                               file_names,
                                               input$selected_node) %>%
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
            
            
            
            plot_damage <- function(x){
                ggplot(x, aes(Position, 
                              Frequency,
                              colour = Mismatch,
                              group = Mismatch)) +
                    geom_line() + 
                    xlab("Position (bp)") +
                    ylab("Alignments (n)") +
                    facet_wrap(File ~ Strand, scales = "free_x", ncol = 4)  + 
                    scale_colour_manual(values = mismatch_colours) +
                    theme_minimal()
            }
            
            
            if (any(selected_filter %in% "default")) {
                damage_plot <- plot_damage(damage_data_comparison %>% 
                                               filter(Mode == "default"))
            } else {
                damage_plot <- plot_damage(damage_dat_comparison)
            }
            
            damage_plot
            
        })
    }
)
