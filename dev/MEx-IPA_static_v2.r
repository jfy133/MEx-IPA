## Pre script stuff

library(tidyverse)
library(data.table)
library(plotly)
library(patchwork)

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
    print(remove_string)
    
    ## Note static requires second str_remove, but not in shiny!
    in_dat %>% 
      mutate(File = str_remove(File, remove_string)) %>%
      filter(File %in% str_remove(!! sample, remove_string),
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

## Load data and get input selection options
select_input <- "~/Documents/Scripts/shiny_web_apps/MEx-IPA/dev/test_data/output_dairymicrobes_archive/"
input_dir <- paste0(select_input,"/")

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
filter_names <- list(all = c("default", "ancient"), 
                     default = c("default"), 
                     ancient = c("ancient"))

## User options selection ######################################################
selected_node <- node_names[37]
selected_file <- file_names[1]
selected_filter <- filter_names[1]  %>% unlist %>% unname
selected_interactive <- T
remove_string <- "_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.rma6"
################################################################################

## Data processing
damage_data <- filter_module_files(damageMismatch, 
                                   selected_filter, 
                                   remove_string,
                                   selected_file,
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
  
edit_data <- filter_module_files(editDistance, 
                                 selected_filter, 
                                 remove_string,
                                 selected_file, 
                                 selected_node) %>%
  gather(Edit_Distance, Alignment_Count, 6:ncol(.)) %>%
  mutate(Edit_Distance = factor(Edit_Distance, levels = c(0:10, "higher")))

percentidentity_data <- filter_module_files(percentIdentity, 
                                            selected_filter, 
                                            remove_string,
                                            selected_file, 
                                            selected_node) %>%
  gather(Percent_Identity, Alignment_Count, 6:ncol(.)) %>%
  mutate(Percent_Identity = factor(Percent_Identity, levels = seq(80, 100, 5)))

length_data <- filter_module_files(readLengthDist, 
                                   selected_filter, 
                                   remove_string, 
                                   selected_file, 
                                   selected_node) %>%
  gather(Length_Bin, Alignment_Count, 6:ncol(.)) %>%
  mutate(Length_Bin = str_replace_all(Length_Bin, "bp", ""),
         Length_Bin = as.numeric(Length_Bin))

positionscovered_data <- filter_module_files(positionsCovered,
                                       selected_filter, 
                                       remove_string,
                                       selected_file, 
                                       selected_node) %>%
  select(-contains("Average"), -contains("_"), -Reference) %>%
  gather(Breadth, Percentage, 6:ncol(.)) %>%
  mutate(Breadth = str_remove_all(Breadth, "percCoveredHigher") %>% 
           as.numeric)

coveragehist_data <- filter_module_files(coverageHist, 
                                    selected_filter, 
                                    remove_string,
                                    selected_file, 
                                    selected_node) %>%
  gather(Fold_Coverage, Base_Pairs, 7:ncol(.)) %>%
  mutate(Fold_Coverage = as_factor(Fold_Coverage)) %>%
  filter(Fold_Coverage != "0")

filterstats_data <- filterStats %>%
  select(-Module) %>%
  left_join(readLengthStat %>% select(-Module), 
            by = c("Index", "Mode", "File", "Node")) %>%
  left_join(alignmentDist %>% select(-Module, -Reference), 
            by = c("Index", "Mode", "File", "Node")) %>% 
  left_join(positionsCovered %>% select(-Module, -Reference, -contains("perc")), 
            by = c("Index", "Mode", "File", "Node")) %>%
  filter_module_files(selected_filter, remove_string, selected_file, selected_node) %>%
  mutate(Mean = round(Mean, digits = 2),
         GeometricMean = round(Mean, digits = 2),
         StandardDev = round(StandardDev, digits = 2)) %>%
  gather(Information, Value, 5:ncol(.)) %>%
  select(-Index) %>% 
  mutate(Information = map(Information, function(x) 
    names(filterstats_info[filterstats_info == x])) %>% unlist) %>%
  mutate(Information = factor(Information, 
                              levels = rev(names(filterstats_info)))) %>%
  arrange(Information)

## Plotting

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


if (any(selected_filter %in% "default")) {
  damage_plot <- plot_damage(damage_data %>% filter(Mode == "default"))
} else {
  damage_plot <- plot_damage(damage_data)
}

length_plot <- plot_col(length_data, 
                        "Length_Bin", 
                        "Alignment_Count",
                        "Read Length Bins (bp)",
                        "Alignments (n)")

edit_plot <- plot_col(edit_data,
                      "Edit_Distance",
                      "Alignment_Count",
                      "Edit Distance",
                      "Alignments (n)")

percentidentity_plot <- plot_col(percentidentity_data,
                                 "Percent_Identity",
                                 "Alignment_Count",
                                 "Sequence Identity (%)",
                                 "Alignments (n)")

positionscovered_plot <- plot_col(positionscovered_data,
                                  "Breadth",
                                  "Percentage",
                                  "Fold Coverage (X)",
                                  "Percentage of Reference (%)")

coveragehist_plot <- plot_col(coveragehist_data,
                              "Fold_Coverage",
                              "Base_Pairs",
                              "Fold Coverage (X)",
                              "Base Pairs (n)")

filterstats_plot <- ggplot(filterstats_data, 
                           aes(y = Information,
                               x = Mode, 
                               label = Value)) +
  geom_tile(colour = "darkgrey", fill = NA, size = 0.3) +
  geom_text() +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank())

## Final plot display

if (!isTRUE(selected_interactive)) {
  damage_plot +
    edit_plot +
    percentidentity_plot +    
    length_plot +
    positionscovered_plot +
    coveragehist_plot
  filterstats_plot
  
} else if (isTRUE(selected_interactive)) {
  ggplotly(damage_plot)
  ggplotly(edit_plot)
  ggplotly(length_plot)
  ggplotly(percentidentity_plot)
  filterstats_plot
}



  




