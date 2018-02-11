library(tidyverse)

zip_list <- unzip("../Documents/Scripts/shiny_web_apps/test_data/90.zip")

file_prefix <- paste(str_split(zip_list[1], "/")[[1]][1:2], collapse="/")

zip_default_runsum <- read_tsv(zip_list[grepl("*/default/RunSummary.txt", zip_list)])%>%
  filter(Node != "Total_Count")
list_sample <- colnames(zip_default_runsum)[2:ncol(zip_default_runsum)] %>% sort
first_taxon <- colnames(zip_default_runsum)[2]
list_taxon <- select(zip_default_runsum, Node) %>%
  distinct()


#dir_damage <- paste(input_dir, "/default/damageMismatch", sep="")

files_damage <- zip_list[grepl("*damageMismatch.txt", zip_list)]

data_damage <- data_frame(filename = files_damage) %>%
  mutate(file_contents = map(filename, ~ read_tsv(.)))

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
#dir_edit <- paste(input_dir, "/default/editDistance", sep="


data_extractor <- function(input, search_key){
  files_out <- zip_list[grepl(search_key, input)]
  data_out <- data_frame(filename = files_out) %>%
    mutate(file_contents = map(filename, ~ read_tsv(.)))
  data_out <- unnest(data_out)
  return(data_out)
}

data_damage <- data_extractor(zip_list, "*_damageMismatch.txt")
data_edit <- data_extractor(zip_list, "*_editDistance.txt")
data_lngt <- data_extractor(zip_list, "*_readLengthDist.txt")



sample <- "Zape5.rma6"
taxon <- "Zea_mays"

final_damage <- data_damage %>%
  filter(grepl(sample, filename)) %>%
  filter(Node == taxon)

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
final_edit <- data_edit %>%
  filter(grepl(sample, filename)) %>%
  filter(Taxon == taxon)

final_edit <- gather(final_edit, Distance, Reads, 3:ncol(final_edit)) %>%
  rowwise() %>%
  mutate(dataset = unlist(str_split(filename, "/"))[3])

edit_reads <- final_edit %>%
  summarise(sum = sum(Reads)) %>%
  sum

distance_levels <- c("0", "1", "2", "3", "4", "5", "6",
                     "7", "8", "9", "10", "higher")

### Read Length
#### Sample and taxon filtering
final_lngt <- data_lngt %>%
  filter(grepl(sample, filename)) %>%
  filter(Taxon == taxon) %>%
  gather(length, reads, 3:ncol(.)) %>%
  mutate(length = str_replace(length, "bp", "")) %>%
  mutate(length = as.numeric(length)) %>%
  rowwise() %>%
  mutate(dataset = unlist(str_split(filename, "/"))[3])

lngt_reads <- final_lngt %>% filter(dataset == "default") %>%
  summarise(sum = sum(reads)) %>%
  sum

## Plot damage
damage_plot <- ggplot(final_damage, aes(x=factor(Position), y=Frequency,
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
                           "reads considered C to T (Red), G to A (Blue)"))
## Plot edit distnace
ggplot() +
  geom_bar(data=filter(final_edit, dataset == "default"), aes(factor(Distance, levels = distance_levels), Reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
  geom_bar(data=filter(final_edit, dataset == "ancient"), aes(factor(Distance, levels = distance_levels), Reads, colour=dataset, fill=dataset), stat="identity", alpha=0.5) +
  xlab("edit distance") +
  ylab("reads") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face="bold")) +
  ggtitle("Edit Distance", subtitle = paste(edit_reads,
                                            "reads considered"))

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
  ggtitle("Read Distribution", subtitle = paste(lngt_reads$sum,
                                                "reads considered"))
## Plot info box
info_text <- paste("Sample:",
                   sample,
                   "\n\n Displayed taxon:",
                   taxon,
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

