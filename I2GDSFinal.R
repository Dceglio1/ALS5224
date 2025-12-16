#Initializing tidyverse, setting working directory, inputting needed files
library(tidyverse)

setwd("C:/Users/Jvlan/Downloads")

#The drugs object comes from the deeparg database, where they assigned each specific antibiotic to a broader class
#The CARD_aro object is from the CARD database, and shows which gene confers resistance to each specific antibiotic
CARD_aro <- read.csv("CARD4.0.1_aro_cat.csv")
Args <- read.delim("ChezLiz_merged_args_new.txt")
metadata <- read.csv("ChezLizMetadata.csv")
drugs <- read.csv("Antibiotic to Drug.csv")

#To make the graph look better, we assign the genes in the CARD database to a broader drug class using the deeparg conversion
for (i in 1:nrow(CARD_aro)) {
  for (j in 1:nrow(drugs)) {
    if(grepl(drugs[j,1], CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- drugs[j,2]}
  }
}

#If a gene confers resistance to multiple drugs, we assign it "multidrug"
for (i in 1:nrow(CARD_aro)) {
  if(grepl(";", CARD_aro[i, 4])) {CARD_aro[i, "Drug"] <- "multidrug"}
}

#If it is still unassigned, it gets classified as other
CARD_aro$Drug[is.na(CARD_aro$Drug)] <- "Other"

#Merging the DIAMOND data and the CARD metadata about each gene
Args <- merge(Args, CARD_aro, by.x = "protein_accession", by.y = "Protein.Accession")

#Renaming the samples to only include the same name
Args$sample = substr(Args$sample, 1, nchar(Args$sample)-10)

#Merging the DIAMOND data and sample metadata
Args <- merge(Args, metadata, by.x = "sample", by.y = "sample_original")

#Just incase the rpob normalization messed up and left this as an NA
#Something went wrong in ARC if this is an NA!
Args$rpob_Normalization[Args$rpob_Normalization==""] <- 0

#Assigning data its appropriate type
Args$Drug.Class <- as.factor(Args$Drug.Class)
Args$Fraction <- as.factor(Args$Fraction)
Args$rpob_Normalization <- as.numeric(Args$rpob_Normalization)
Args$count <- as.numeric(Args$count)
Args$Drug <- as.factor(Args$Drug)

#Creating df where bubble plot will be generated from
#Fraction = Sampling location
#Only separating by fraction to see how different sampling locations are different
bubble <- Args %>%
  group_by(Drug, Fraction) %>%
  summarize(sum_rpob_normalized_count = mean(rpob_Normalization))

#The "CONTROL" sample type occurs every 5th row in bubble
controlRow <- 5

#Creating new column in bubble to put log differences
bubble$logDiff <- 0

#This loop goes through the df and calculates the log difference in 
#rpob normalized ARGs compared to the control for each drug class
for (i in 1:13) { #Number of drug classes
  for (k in 1:7) { #Number of sample types
    count = k+((i*7)-7)
    bubble[count, 4] <- log10(bubble[count, 3]/bubble[controlRow, 3])
  }
  controlRow <- controlRow +7
}

#Find absolute value for plotting
bubble$absLogDiff <- abs(bubble$logDiff)

#Differentiate between increased ARGs (logDiff > 0) and decreased (logDiff < 0)
bubble$Pos <- ifelse(bubble$logDiff >= 0, "Increase", "Decrease")

#Set order of Fraction for graph
bubble$Fraction <- factor(bubble$Fraction, levels = c("INF", "EFF", "BOIL", "30M", "100M", "Blank", "CONTROL"))

#Using ggplot to plot the graph
ggplot(subset(bubble, Fraction %in% c("100M", "30M", "BOIL", "EFF", "INF")), aes(x = Fraction, y = Drug, size = absLogDiff)) +
  geom_point(alpha = 0.8, aes(color = Pos)) +
  theme_bw() + 
  scale_size(range = c(1, 20)) +
  xlab("Sampling Location") +
  labs(size = "Log Difference to Control") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(color = "Change") 
  





