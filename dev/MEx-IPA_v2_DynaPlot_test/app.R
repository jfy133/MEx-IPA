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

## Load data specific information for user options
node_names <- default_runsummary %>% pull(Node) 
file_names <- default_runsummary %>% colnames %>% .[. != "Node"]
filter_names <- list("all", "default", "ancient")

remove_string <- ""
selected_node <- "Lactobacillus_ruminis"
selected_filter <- "default"


damage_data_comparison <- filter_module_files(damageMismatch, 
                                              selected_filter, 
                                              "",
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


# Define UI for application that draws a histogram
ui <- fluidPage(
    numericInput("number", label = NULL, value = 1, step = 1, min = 1),
    uiOutput("plots")
)


# Define server logic required to draw a histogram
server <- function(input, output) {

    # plotInput <- reactive({
    #     n_plot <- input$number
    #     total_data <- lapply(1:n_plot, function(i){rnorm(500)})
    #     print(total_data)
    #     return(list("n_plot" = n_plot, "total_data" = total_data))
    #     print(list("n_plot" = n_plot, "total_data" = total_data))
    # })
    
    plotInput <- reactive({
        n_plot <- file_names
        total_data <- map(n_plot, function(x) {
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
            
        })
        names(total_data) <- n_plot
        return(list("n_plot" = n_plot, "total_data" = total_data))
    }) 
    
    
    ##### Create divs######
    output$plots <- renderUI({
        plot_output_list <- lapply(plotInput()$n_plot, function(i) {
            plotname <- i
            plotOutput(plotname)
        })   
        #print(plot_output_list)
        do.call(tagList, plot_output_list)
    })
    
    observe({
        lapply(plotInput()$n_plot, function(i){
            output[[i]] <- renderPlot({
                plot_damage(plotInput()$total_data[[i]])
            })
        })
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
