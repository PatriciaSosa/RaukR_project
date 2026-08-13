library(data.table)
library(tibble)
library(tidyverse)
library(dplyr)
library(tidyr)
library(corrplot)
library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
# BiocManager::install("factoextra")
library(factoextra)
library(ggplot2)
# BiocManager::install("umap")
library(umap)
library(readxl)
library(ggrepel)
library(ggpmisc)

# color palette for the QC data ---
# blank <- "#0E131F"
# pool <- "#077187"
# sample <- "#FCA17D"

####### Serology experiment for Uganda encephalitis samples ########
# Loading MFI values, antigens, sample position, and metadata from CSV files
#data <- read.csv("PBA012_plate1_ENC_median.csv", header = TRUE) # controls plate 1 + cases plate 1 + cases plate 2 + cases plate 3 (first run)
data2 <- read.csv("PBA012_plate2_ENC_median.csv", header = TRUE) # controls plate 2

data <- read.csv("PBA012_RERUN_ENC_median.csv", header = TRUE) # controls plate 1 + cases plate 1 + cases plate 2 + cases plate 3

positions <- read_xlsx("384_plate_layout_traceability_v2.xlsx")
positions2 <- read_xlsx("Coari_Itacoatiara_Plate3B_96w_to_384w_traceability.xlsx")

antigens <- read_xlsx("antigen_dilution_ENC_run2.xlsx")

# removing columns and rows from the antigens data frame that are not needed
antigens <- cbind(antigens$Name,
                  antigens$bead_id,
                  antigens$Antigen_Group,
                  antigens$Pathogen,
                  antigens$virus_genus,
                  antigens$protein)

antigens <-as.data.frame(antigens)
colnames(antigens) <- c( "antigen", "bead_id", "antigen_group","Pathogen",
                         "virus_genus", "protein")


##### Converting wide format to long format with the RAW MFI VALUES ######!!!!!!!
data_long <- data %>% pivot_longer(cols = -c(Location, dest_well, Sample, Total.Events),
                                       names_to = "antigen",
                                       values_to = "MFI")

data_long2 <- data2 %>% pivot_longer(cols = -c(Location, dest_well, Sample, Total.Events),
                                   names_to = "antigen",
                                   values_to = "MFI")

# creating a new column in data_long separating antigen numbers to a new column
data_long <- data_long %>%
  separate(antigen, into = c("antigen", "bead_id"), sep = "Analyte.")

data_long2 <- data_long2 %>%
  separate(antigen, into = c("antigen", "bead_id"), sep = "Analyte.")

# removing antigen column from data_long1
data_long <- data_long %>% select(-antigen)
data_long2 <- data_long2 %>% select(-antigen)

# merging data_long1 with antigens data frame to get antigen information
merged_data <- merge(data_long, antigens, by = c("bead_id"))
merged_data2 <- merge(data_long2, antigens, by = c("bead_id"))

# merging merged_data with positions data frame to get sample metadata
preliminar_data <- merge(merged_data, positions, by.x = "dest_well", by.y = "well_384")
preliminar_data2 <- merge(merged_data2, positions2, by.x = "dest_well", by.y = "well_384")

# filtering preliminar_data by cohort to remove unwanted cohorts
preliminar_data2 <- preliminar_data2 %>% 
  filter(selma_quadrant == "Q4")

preliminar_data <- cbind(preliminar_data$dest_well,
                         preliminar_data$bead_id,
                         preliminar_data$Sample,
                         preliminar_data$MFI,
                         preliminar_data$antigen,
                         preliminar_data$antigen_group,
                         preliminar_data$Pathogen,
                         preliminar_data$virus_genus,
                         preliminar_data$protein,
                         preliminar_data$sample_id,
                         preliminar_data$group)

preliminar_data <-as.data.frame(preliminar_data)
colnames(preliminar_data) <- c("dest_well", "bead_id", "Sample", "MFI", "antigen", "antigen_group", "Pathogen",
                         "virus_genus", "protein", "sample_id", "group")

preliminar_data2 <- cbind(preliminar_data2$dest_well,
                          preliminar_data2$bead_id,
                          preliminar_data2$Sample,
                          preliminar_data2$MFI,
                          preliminar_data2$antigen,
                          preliminar_data2$antigen_group,
                          preliminar_data2$Pathogen,
                          preliminar_data2$virus_genus,
                          preliminar_data2$protein,
                          preliminar_data2$sample_id,
                          preliminar_data2$group)

preliminar_data2 <-as.data.frame(preliminar_data2)
colnames(preliminar_data2) <- c("dest_well", "bead_id", "Sample", "MFI", "antigen", "antigen_group", "Pathogen",
                               "virus_genus", "protein", "sample_id", "group")

# creating a plate column in preliminar_data and preliminar_data2 to indicate which plate the data is from
preliminar_data$plate <- "plate1"
preliminar_data2$plate <- "plate2"

# combining preliminar_data and preliminar_data2 into a single data frame
ENC02_data <- rbind(preliminar_data, preliminar_data2)
ENC02_data$MFI <- as.numeric(ENC02_data$MFI)

preliminar_data$MFI <- as.numeric(preliminar_data$MFI)
preliminar_data2$MFI <- as.numeric(preliminar_data2$MFI)

preliminar_data$antigen <- reorder(preliminar_data$antigen,
                                   preliminar_data$MFI, median)

preliminar_data2$antigen <- reorder(preliminar_data2$antigen,
                                   preliminar_data2$MFI, median)

ENC02_data$antigen <- reorder(ENC02_data$antigen,
                              ENC02_data$MFI,
                              median)

######### checking data quality ############
ggplot(preliminar_data %>% filter(virus_genus == "NA"), aes(x = antigen, y = MFI)) +
  geom_jitter(aes(color = group), width = 0.2) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by QC beads (plate 1)",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  facet_wrap(~plate, scales = "free")

ggplot(preliminar_data %>% filter(virus_genus == "NA"), aes(x = antigen, y = MFI)) +
  geom_jitter(aes(color = group), width = 0.2) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by QC beads (plate 1)",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  facet_wrap(~plate, scales = "free") +
  ylim(0, 2500)

ggplot(preliminar_data %>% filter(group == "blank"), aes(x = antigen, y = MFI)) +
  geom_jitter(aes(color = Pathogen),width = 0.2) +
  theme_bw() +
  labs(title = "MFI values in blank samples",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  facet_wrap(~plate, scales = "free")


ggplot(preliminar_data %>% filter(group == "blank"), aes(x = antigen, y = MFI)) +
  geom_violin() +
  geom_jitter(aes(color = Pathogen),width = 0.2) +
  theme_bw() +
  labs(title = "MFI values in blank samples",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  facet_wrap(~plate, scales = "free") +
  ylim(0, 100)

# both plates 1 and 2
ggplot(ENC02_data %>% filter(virus_genus == "NA"), aes(x = antigen, y = MFI)) +
  geom_jitter(aes(color = group), width = 0.2) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by QC beads",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  facet_wrap(~plate, scales = "free")

ggplot(ENC02_data, aes(x = antigen, y = MFI, fill = antigen_group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by Antigen",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15))

ggplot(ENC02_data, aes(x = antigen, y = MFI, fill = group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by Antigen",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15))

ggplot(ENC02_data, aes(x = antigen, y = MFI, fill = antigen_group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by Antigen",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15))

ggplot(ENC02_data, aes(x = antigen, y = MFI, fill = antigen_group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by Antigen",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  ylim(0, 5000)

ggplot(ENC02_data, aes(x = antigen, y = MFI, fill = antigen_group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("case" = "#B6174B", "control" = "#FCA17D", "blank" = "#0E131F", "pool" = "#077187")) +
  theme_bw() +
  labs(title = "MFI values by Antigen",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 14), 
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  ylim(0, 5000)

ggplot(ENC02_data, aes(x = antigen, y = MFI, fill = antigen_group)) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  labs(title = "MFI values by Antigen",
       x = "Antigens",
       y = "Raw MFI Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 10), 
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  facet_wrap(~group, scales = "free")


