#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
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

clean_damage <- function(x) {
    filter_module_files(damageMismatch, 
                        selected_filter, 
                        remove_string,
                        x,
                        selected_node) %>%
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

clean_length <- function(x) {
    filter_module_files(readLengthDist, 
                        selected_filter, 
                        remove_string,
                        x,
                        selected_node) %>%
        gather(Length_Bin, Alignment_Count, 6:ncol(.)) %>%
        mutate(Length_Bin = str_replace_all(Length_Bin, "bp", ""),
               Length_Bin = as.numeric(Length_Bin))
}

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

## default aesthetics

mode_colours <- c(default = "#1b9e77", ancient = "#7570b3") 

mismatch_colours <- c(`C>T` = "#e41a1c", 
                      `G>A` = "#377eb8",
                      `D>V(11Substitution)` = "grey", 
                      `H>B(11Substitution)` = "grey")

damage_xaxis <- c("11" = "-10", "12" = "-9", "13" = "-8", "14" = "-7",
                  "15" = "-6", "16" = "-5", "17" = "-4", "18" = "-3",
                  "19" = "-2", "20" = "-1")

## Load data

input_dir <- "/home/fellows/Documents/Scripts/shiny_web_apps/MEx-IPA/dev/test_data/output_dairymicrobes_archive/"
default_runsummary <- load_runsummary(input_dir, "default/")
damageMismatch <- load_module_files(input_dir, "_damageMismatch.txt")
readLengthDist <- load_module_files(input_dir, "_readLengthDist.txt")

## Load data specific information for user options
node_names <- default_runsummary %>% pull(Node) 
file_names <- default_runsummary %>% colnames %>% .[. != "Node"]
filter_names <- list("all", "default", "ancient")

remove_string <- ""
selected_node <- "Lactobacillus_ruminis_ATCC_27782"
selected_filter <- "default"
interactive <- F
statistic <- "length" # or "length"

damage_data_comparison <- filter_module_files(damageMismatch, 
                                              selected_filter, 
                                              remove_string,
                                              file_names,
                                              selected_node) %>%
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

length_data_comparison <- filter_module_files(readLengthDist, 
                                              selected_filter, 
                                              remove_string,
                                              file_names,
                                              selected_node) %>%
    gather(Length_Bin, Alignment_Count, 6:ncol(.)) %>%
    mutate(Length_Bin = str_replace_all(Length_Bin, "bp", ""),
           Length_Bin = as.numeric(Length_Bin))


# Define UI for application that draws a histogram
ui <- fluidPage(
    uiOutput("plots")
)


# Define server logic required to draw a histogram
server <- function(input, output) {
    
    ## 1) Make a list of all the data for each plot
    ## 2) Make placeholders for each plot's data listed in the list
    ## 3) For each entry in list, render the plot (observe executes immediately
    ## when it detects a change)
    ## Notes: as using renderUI, sends all plots to UI at once.
    
    plotInput <- reactive({
        n_plot <- file_names
        
        if (statistic == "damage") {
            total_data <- map(n_plot, clean_damage)
        } else if (statistic == "length") {
            total_data <- map(n_plot, clean_length)
        }
        
        names(total_data) <- n_plot
        
        
        ## Removes entries with empty tibbles, n_plot above will be for every 
        ## file, even if no data
        total_data <- total_data[map(total_data, nrow)>0]
        
        ## reset n_plot
        n_plot <- names(total_data)
        
        return(list("n_plot" = n_plot, "total_data" = total_data))
    })
    
    
    ##### Create divs######
    output$plots <- renderUI({
        
        ## make all the plot(ly) objects and place in a list
        if (interactive) {
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
        if (interactive) {
            lapply(plotInput()$n_plot, function(i){
                output[[i]] <- renderPlotly({
                    if (statistic == "damage") {
                        plot_damage(plotInput()$total_data[[i]])
                    } else if (statistic == "length") {
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
                    if (statistic == "damage") {
                        plot_damage(plotInput()$total_data[[i]])
                    } else if (statistic == "length") {
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
}

# Run the application 
shinyApp(ui = ui, server = server)
