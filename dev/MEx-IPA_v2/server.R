#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(data.table)
library(plotly)
library(patchwork)
library(shinyFiles)

## Functions 

load_runsummary <- function(in_dir, result_path) {
  fread(paste0(in_dir, result_path, "RunSummary.txt")) %>% 
    as_tibble
}


load_module_files <- function(in_dir, file_ext) {
  list.files(input_dir, 
             pattern = file_ext, 
             recursive = T,
             full.names = T) %>% 
    enframe("Index", "File") %>%
    mutate(data = purrr::map(File, function(x) fread(x, header = T) %>% as_tibble),
           File = str_remove_all(File,
                                 paste(c(in_dir, 
                                         file_ext, "^/"), 
                                       collapse = "|"))) %>% 
    separate(File, 
             c("empty", "Mode", "Module", "File"), 
             sep =  "/") %>% 
    select(-empty) %>%
    unnest %>%
    rename_at(vars(matches("Taxon")), 
              function(x) str_replace_all(x, "Taxon", "Node")) %>%
    mutate(Mode = factor(Mode, levels = c("default", "ancient")))
}

filter_module_files <- function(in_dat, filt, File, taxon) {
  sample <- enquo(File)
  taxon <- enquo(taxon)
  in_dat %>% 
    filter(File == !! sample,
           Node == !! taxon,
           Mode %in% filt)
}

plot_damage <- function(x){
  ggplot(x, aes(Position, 
                Frequency,
                colour = Mismatch,
                group = Mismatch)) +
    geom_line() + 
    labs(title = "Misincorporation (Damage) Plot") +
    xlab("Position (bp)") +
    ylab("Alignments (n)") +
    facet_wrap(~ Strand, scales = "free_x")  + 
    scale_colour_manual(values = mismatch_colours) +
    theme_minimal()
}

## Aesthetics

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

# Server Logic
shinyServer(function(input, output) {
   
  shinyDirChoose(input, 
                 "input_dir", 
                 updateFreq = 1, 
                 roots = c(home = '~'),
                 defaultPath = ".", 
                 defaultRoot = "~")
  
  input_dir <- reactive(input$input_dir)
  
  output$input_dir <- renderText({ 
    paste("You have selected this", input$input_dir)
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 home <- normalizePath("~")
                 datapath <<-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
               })
  
})
