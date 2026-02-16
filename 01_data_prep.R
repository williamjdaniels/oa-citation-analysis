# ==============================================================================
# Script 1: Data Preparation
# Citation Analysis of Open Access Articles
# 
# Input:  Raw CSV data files
# Output: min_data.csv (processed data ready for analysis)
# ==============================================================================

# Load required packages
library(dplyr)
library(tidyverse)
library(vroom)

# Functions
count_commas <- function(text) {
  return(lengths(regmatches(text, gregexpr(",", text))))
}

# Load and prepare data
data <- vroom::vroom("All Journals and Articles with Stratified OA.xlsx - Sheet1.csv")
data$`Open Access Stratified` <- ifelse(data$`Open Access Stratified` == "Hybrid", 
                                        "Gold", 
                                        data$`Open Access Stratified`)

# Load AIS data
AIS_dat <- vroom::vroom("~/Downloads/AUCompBio-OA_2021_AU-cd94c12/AIS2023_dat - Sheet1.csv")

# Select relevant variables
min_data <- dplyr::select(data, 
                          "Journal Title", 
                          "Publication Year", 
                          "Open Access Stratified", 
                          "Times Cited, All Databases",
                          "Authors",
                          "Document Type")

Articles_Per_Year <- min_data %>% group_by(`Publication Year`) %>% summarize(n())
message(paste0("Dropping ", 
               sum(Articles_Per_Year$`n()`[7:9]), 
               " Articles of ", 
               sum(Articles_Per_Year$`n()`[1:9]), " Remaining articles: ", sum(Articles_Per_Year$`n()`[1:6])))

# Remove recent years with incomplete citation windows
min_data <- subset(min_data, !`Publication Year` %in% c(2020:2022))
doc_type <- min_data %>%
  group_by(`Document Type`) %>%
  summarize(num_type = n()) %>%
  arrange(num_type)
tail_doc <- tail(doc_type, 6)
min_data <- subset(min_data, `Document Type` %in% tail_doc$`Document Type`)

# Calculate author counts
min_data$comma_count <- count_commas(min_data$`Authors`)
min_data$author_count <- min_data$comma_count + 1

# Load and merge AIS data
min_data <- left_join(min_data, dplyr::select(AIS_dat, AIS, `Journal Title`))

# Handle missing AIS values - antequated due to NO MISSING AIS SCORES!
# missing_ais <- sum(is.na(min_data$AIS))
# if(missing_ais > 0) {
#   min_data$AIS[is.na(min_data$AIS)] <- mean(min_data$AIS, na.rm = TRUE)
# }

# Scale continuous variables (preserving parameters for back-transformation)
author_mean <- mean(min_data$author_count, na.rm = TRUE)
author_sd <- sd(min_data$author_count, na.rm = TRUE)
min_data$auth_count_scaled <- (min_data$author_count - author_mean) / author_sd

AIS_mean <- mean(min_data$AIS, na.rm = TRUE)
AIS_sd <- sd(min_data$AIS, na.rm = TRUE)
min_data$AIS_scaled <- (min_data$AIS - AIS_mean) / AIS_sd

year_mean <- mean(min_data$`Publication Year`, na.rm = TRUE)
year_sd <- sd(min_data$`Publication Year`, na.rm = TRUE)
min_data$pub_year_scaled <- (min_data$`Publication Year` - year_mean) / year_sd

# Factor variables with correct order
min_data$`Open Access Stratified` <- factor(min_data$`Open Access Stratified`, 
                                            levels = c("Bronze", "Closed Access", "Green", "Gold"))
min_data$`Journal Title` <- factor(min_data$`Journal Title`)

# Create squared term for potential nonlinear effects
min_data$pub_year_scaled_sq <- min_data$pub_year_scaled^2

# Save scaling parameters for use in analysis script
scaling_params <- data.frame(
  variable = c("author_count", "AIS", "Publication Year"),
  mean = c(author_mean, AIS_mean, year_mean),
  sd = c(author_sd, AIS_sd, year_sd)
)
write.csv(scaling_params, "scaling_params.csv", row.names = FALSE)

## Write data for publication
data.table::fwrite(min_data, "min_data.csv")

cat("\n=== Data Preparation Complete ===\n")
cat(sprintf("Total articles in processed dataset: %d\n", nrow(min_data)))
cat(sprintf("Number of journals: %d\n", n_distinct(min_data$`Journal Title`)))
cat(sprintf("Year range: %d-%d\n", min(min_data$`Publication Year`), max(min_data$`Publication Year`)))
cat(sprintf("Mean authors per article: %.1f (SD = %.1f)\n", author_mean, author_sd))
cat("Output saved: min_data.csv, scaling_params.csv\n")