# Clear Workspace
rm(list = ls())

# Necessary libraries
library(caret)
library(plyr)
library(Metrics)
library(mice)
library(ggplot2)
library(RColorBrewer)

# Set seed for reproducibility and also set working directory
set.seed(1)
setwd("C:/Users/User/Documents/R/BROADPEAKSR/")


score_threshold = 600

# Load the  data
raw_data <- read.csv(file = "h3k4me3_rep1_peaks.bed", sep = "\t", col.names = c("chrom", "start", "end", "score"))

islands_above_score <- matrix(ncol=2, nrow=score_threshold)

for (i in 1:score_threshold ) {
        islands_above_score[i,1] = i
        islands_above_score[i,2] = sum(raw_data$score >= i)
}

islands_above_score = data.frame(islands_above_score)
colnames(islands_above_score) = c("score_thres", "number_of_islands")
qplot(score_thres, number_of_islands, data = islands_above_score)
