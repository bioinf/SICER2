# Clear Workspace
rm(list = ls())

# Necessary libraries
library(caret)
library(plyr)
library(Metrics)
library(mice)
library(ggplot2)
library(RColorBrewer)
library(rms)

# Set seed for reproducibility and also set working directory
set.seed(1)
setwd("C:/Users/User/Documents/R/BROADPEAKSR/")

# Load the  data
raw_data <- read.csv(file = "e_val.txt", sep = "\t", col.names = c("e_val", "score_thres", "p_val", "number_of_reads"))
inTrain         = createDataPartition(raw_data$score_thres, p = 0.6)[[1]]
training        <- raw_data[inTrain,]
remainder        <- raw_data[-inTrain,]


model <- glm(score_thres~(e_val) + p_val + number_of_reads,  
              data = training)


ggplot(raw_data, aes(x=p_val, y=score_thres, color = log(e_val))) + geom_point(shape=1)+ geom_smooth(method=lm)  
ggplot(raw_data, aes(x=log(e_val), y=score_thres, color = p_val)) + geom_point(shape=1)+ geom_smooth(method=lm)
