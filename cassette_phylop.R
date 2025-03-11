library(tidyverse)
library("patchwork")
library("ggpubr")
library("rtracklayer")
library("GenomicRanges")

# import .bw (conservation score) file:
bw <- import.bw("ce11.phyloP26way.bw")

# Extract start and end coordinates of each cassette exon event:

extracted_coords <- lapply(rownames(cassette_sums_heatmap), function(id) {
  event_id <- sub("^[^ ]+ ", "", id)
  parts <- unlist(strsplit(event_id, "_"))
  chr <- as.character(parts[1])
  region_start <- as.numeric(parts[3])
  cassette_start <- as.numeric(parts[4])
  cassette_end <- as.numeric(parts[5])
  region_end <- as.numeric(parts[6])
  return(c(chr, region_start, cassette_start, cassette_end, region_end))
})

extracted_coords <- do.call(rbind, extracted_coords)
rownames(extracted_coords) <- sub("^[^ ]+ ", "", rownames(cassette_sums_heatmap))
#write.csv(extracted_coords, file = "extracted_coords_cassette_events.csv")   ### move this file to HPCC and run average_phylop_score_for_cassette_exons.txt


# calculate mean phylop scores between extracted coordinates ranges

# Initialize output data frame
phyloP_scores_combined <- data.frame(AS_event_ID = character(), Chr = character(), region_start = integer(), cassette_start = integer(), cassette_end = integer(), region_end = integer(), Average_Score = numeric(), stringsAsFactors = FALSE)

# Process each row of extracted_coords
for (i in seq_len(nrow(extracted_coords))) {
  AS_event_ID <- rownames(extracted_coords)[i]
  chr <- paste0("chr", extracted_coords[i, 1])
  region_start <- as.integer(extracted_coords[i, 2])
  cassette_start <- as.integer(extracted_coords[i, 3])
  cassette_end <- as.integer(extracted_coords[i, 4])
  region_end <- as.integer(extracted_coords[i, 5])
  
  # Check if start and end are valid
  if (!is.na(region_start) && !is.na(region_end) && region_start < region_end &&
      !is.na(cassette_start) && !is.na(cassette_end) && cassette_start < cassette_end) {
    # Extract phyloP scores from the BigWig file
    Total_Region_Score <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(region_start, region_end)))$score
    Cassette_Exon_Score <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(cassette_start, cassette_end)))$score
    Flanking_Regions_Score_1 <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(region_start, cassette_start)))$score
    Flanking_Regions_Score_2 <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(cassette_end, region_end)))$score
    Flanking_Regions_Score <- c(Flanking_Regions_Score_1, Flanking_Regions_Score_2)
    
    # Compute average scores
    Total_Region_Score <- if (length(Total_Region_Score) > 0) mean(Total_Region_Score, na.rm = TRUE) else NA
    Cassette_Exon_Score <- if (length(Cassette_Exon_Score) > 0) mean(Cassette_Exon_Score, na.rm = TRUE) else NA
    Flanking_Regions_Score <- if (length(Flanking_Regions_Score) > 0) mean(Flanking_Regions_Score, na.rm = TRUE) else NA
    
  # Debugging output
  cat(sprintf("Processing: AS_event_ID=%s, Chr=%s, Start=%d, End=%d, Total Region Score=%f\n", AS_event_ID, chr, region_start, region_end, Total_Region_Score))
  
    # Append result to output data frame
    phyloP_scores_combined <- rbind(phyloP_scores_combined, data.frame(AS_event_ID, chr, region_start, cassette_start, cassette_end, region_end, Total_Region_Score, Cassette_Exon_Score, Flanking_Regions_Score, stringsAsFactors = FALSE))
  } else {
    cat(sprintf("Skipping invalid coordinates: AS_event_ID=%s (%s, %d, %d)\n", AS_event_ID, chr, region_start, region_end))
  }
}


# Modify and format phyloP_scores_combined into cassette_phylop:

cassette_phylop <- as.data.frame(phyloP_scores_combined)
cassette_phylop$Total_Region_Score <- as.numeric(cassette_phylop$Total_Region_Score)   # convert the score columns to numeric
cassette_phylop$Cassette_Exon_Score <- as.numeric(cassette_phylop$Cassette_Exon_Score)   # convert the score columns to numeric
cassette_phylop$Flanking_Regions_Score <- as.numeric(cassette_phylop$Flanking_Regions_Score)   # convert the score columns to numeric


# Correlation of uniqueness index values with PhyloP score:

cassette_sums_stats <- cassette_sums_heatmap %>%
  rownames_to_column("AS_event_ID") %>%
  mutate(
    Min_Value = pmin(!!!syms(singletypes), na.rm = TRUE),
    Max_Value = pmax(!!!syms(singletypes), na.rm = TRUE)) %>%
  mutate(
    across(c(Min_Value, Max_Value), 
           ~ cut(.x, 
                 breaks = quantile(.x, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                 labels = c("0-20", "21-40", "41-60", "61-80", "81-100"), 
                 include.lowest = TRUE), 
           .names = "{.col}_Percentile"))

cassette_sums_stats <- cassette_sums_stats %>%
  mutate(AS_event_ID = sub("^[^ ]+ ", "", AS_event_ID))
combined_data <- merge(cassette_phylop, cassette_sums_stats, by.y = "AS_event_ID")
combined_data <- combined_data[!is.na(combined_data$Cassette_Exon_Score),]

create_regression_plot <- function(data, x_var, y_var = "Cassette_Exon_Score") {
  plot_data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  
  lm_model <- lm(as.formula(paste(y_var, "~", x_var)), data = plot_data)
  model_summary <- summary(lm_model)
  
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  r_squared <- model_summary$adj.r.squared
  p_value <- coef(model_summary)[2,4]
  
  custom_label <- paste0(
    "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
    "\nAdj. R² = ", round(r_squared, 3),
    "\nP = ", format.pval(p_value, digits = 3))
  
  ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = "dodgerblue", size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", 
             x = quantile(plot_data[[x_var]], 0.7, na.rm = TRUE), 
             y = quantile(plot_data[[y_var]], 0.9, na.rm = TRUE), 
             label = custom_label, size = 5, color = "black", hjust = 0) +
    labs(
      x = paste(gsub("_", " ", x_var), "Uniqueness Value"), 
      y = "Cassette Exon PhyloP Score", 
      title = paste("Linear Regression:", gsub("_", " ", x_var), "Uniqueness Value vs Cassette Exon PhyloP Score")
    ) +
    theme_minimal()
}

max_value_plot <- create_regression_plot(combined_data, "Max_Value")
print(max_value_plot)
min_value_plot <- create_regression_plot(combined_data, "Min_Value")
print(min_value_plot)

# Optional: Combine plots if needed
max_value_plot + min_value_plot

## Absolute value of uniqueness index:
cassette_sums_stats <- cassette_sums_heatmap %>%
  rownames_to_column("AS_event_ID") %>%
  mutate(
    Min_Value_Abs = abs(pmin(!!!syms(singletypes), na.rm = TRUE)),
    Max_Value_Abs = abs(pmax(!!!syms(singletypes), na.rm = TRUE)),
    Highest_Abs = pmax(Min_Value_Abs, Max_Value_Abs)) %>%
  mutate(
    across(c(Highest_Abs), 
           ~ cut(.x, 
                 breaks = quantile(.x, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                 labels = c("0-20", "21-40", "41-60", "61-80", "81-100"), 
                 include.lowest = TRUE), 
           .names = "{.col}_Percentile"))

cassette_sums_stats <- cassette_sums_stats %>%
  mutate(AS_event_ID = sub("^[^ ]+ ", "", AS_event_ID))
combined_data <- merge(cassette_phylop, cassette_sums_stats, by.y = "AS_event_ID")
combined_data <- combined_data[!is.na(combined_data$Cassette_Exon_Score),]

create_regression_plot <- function(data, x_var, y_var = "Cassette_Exon_Score") {
  plot_data <- data %>%
    filter(!is.na(.data[[x_var]]), 
           !is.na(.data[[y_var]]),
           is.finite(.data[[x_var]]),
           is.finite(.data[[y_var]]))
  
  if(nrow(plot_data) < 10) {
    warning("Not enough valid data points for regression")
    return(NULL)
  }
  
  lm_model <- lm(as.formula(paste(y_var, "~", x_var)), data = plot_data)
  model_summary <- summary(lm_model)
  
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  r_squared <- model_summary$adj.r.squared
  p_value <- coef(model_summary)[2,4]
  
  custom_label <- paste0(
    "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
    "\nAdj. R² = ", round(r_squared, 3),
    "\nP = ", format.pval(p_value, digits = 3))
  
  ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = "dodgerblue", size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", 
             x = quantile(plot_data[[x_var]], 0.7, na.rm = TRUE), 
             y = quantile(plot_data[[y_var]], 0.9, na.rm = TRUE), 
             label = custom_label, size = 5, color = "black", hjust = 0) +
    labs(
      x = paste(gsub("_", " ", x_var), "Uniqueness Value"), 
      y = "Cassette Exon PhyloP Score", 
      title = paste("Linear Regression:", gsub("_", " ", x_var), "Uniqueness Value vs Cassette Exon PhyloP Score")) +
    theme_minimal()
}

# Add additional data cleaning before regression
combined_data <- combined_data %>%
  mutate(
    Highest_Abs = as.numeric(as.character(Highest_Abs)),
    Cassette_Exon_Score = as.numeric(as.character(Cassette_Exon_Score))) %>%
  filter(!is.na(Highest_Abs), 
         !is.na(Cassette_Exon_Score),
         is.finite(Highest_Abs),
         is.finite(Cassette_Exon_Score))

highest_uniqueness_plot <- create_regression_plot(combined_data, "Highest_Abs")
print(highest_uniqueness_plot)


# In the last regression, which data points represent cassette exons with high abs. uniqueness values AND high conservation score?

# Define thresholds:
uniqueness_threshold <- quantile(combined_data$Highest_Abs, 0.9, na.rm = TRUE)   ### change threshold here
phylop_threshold <- quantile(combined_data$Cassette_Exon_Score, 0.9, na.rm = TRUE)   ### change threshold here

# Filter and arrange
high_uniqueness_high_phylop <- combined_data %>%
  filter(
    Highest_Abs >= uniqueness_threshold,
    Cassette_Exon_Score >= phylop_threshold) %>%
  arrange(desc(Highest_Abs), desc(Cassette_Exon_Score))

high_uniqueness_high_phylop <- high_uniqueness_high_phylop %>%
  dplyr::select(AS_event_ID, Highest_Abs, Cassette_Exon_Score, everything())

# Optional - add gene name for each AS_event_ID:
high_uniqueness_high_phylop <- high_uniqueness_high_phylop %>%
  rowwise() %>%
  mutate(
    Gene = mega_cassette$Gene[mega_cassette$AS_event_ID == .data$AS_event_ID][1]) %>%
  ungroup()

# Optional - reorder columns to insert Gene between Chr and Region_Start
high_uniqueness_high_phylop <- high_uniqueness_high_phylop[, c(
  "AS_event_ID", 
  "Highest_Abs", 
  "Cassette_Exon_Score", 
  "chr", 
  "Gene",
  "region_start", 
  setdiff(names(high_uniqueness_high_phylop), c("AS_event_ID", "Highest_Abs", "Cassette_Exon_Score", "chr", "region_start", "Gene")))]

# Return the number of events meeting the criteria
print(paste("Number of cassette exons above high uniqueness, high phyloP threshold:", nrow(high_uniqueness_high_phylop)))

# Optional - view the top cassette events:
head(high_uniqueness_high_phylop, 20)
# optional - write.csv of results:
write.csv(head(high_uniqueness_high_phylop, 20), file = "cassette_high_uniqueness_high_phylop.csv")

# Optional - create a scatter plot highlighting these events:
highlight_plot <- ggplot(combined_data, aes(x = Highest_Abs, y = Cassette_Exon_Score)) +
  geom_point(color = "gray", alpha = 0.5) +
  geom_point(data = high_uniqueness_high_phylop, color = "red2", size = 2.8) +
  geom_hline(yintercept = phylop_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = uniqueness_threshold, linetype = "dashed", color = "blue") +
  labs(
    title = "Cassette Exons with High (Abs.) Uniqueness and High PhyloP Score",
    x = "Highest Absolute Uniqueness Value (any cell type)",
    y = "Cassette Exon PhyloP Score") +
  theme_minimal()
print(highlight_plot)


# Do the same linear regression, but alter it for unique neuron types:

regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)

# Create a directory to store plots and initialize empty list to store them in R:
dir.create("uniqueness_regression_plots", showWarnings = FALSE)
plot_list <- list()

for (i in singletypes) {
  uniqueness_values <- as.numeric(cassette_sums_heatmap[[i]])
  regression_data <- data.frame(
    AS_event_ID = rownames(cassette_sums_heatmap),
    uniqueness_values = uniqueness_values,
    stringsAsFactors = FALSE)
  
  regression_data$AS_event_ID <- sub("^[^ ]+ ", "", regression_data$AS_event_ID)
  
  merged_data <- merge(regression_data, cassette_phylop, by = "AS_event_ID", all.x = TRUE)
  regression_data <- merged_data[!is.na(merged_data$Cassette_Exon_Score),]
  regression_data <- regression_data[regression_data$uniqueness_values != 0,]
  regression_data <- regression_data[regression_data$Cassette_Exon_Score != "No Data",]
  regression_data$Total_Region_Score <- as.numeric(regression_data$Total_Region_Score)
  regression_data$Cassette_Exon_Score <- as.numeric(regression_data$Cassette_Exon_Score)
  regression_data$Flanking_Regions_Score <- as.numeric(regression_data$Flanking_Regions_Score)
  
  lm_model <- lm(Cassette_Exon_Score ~ uniqueness_values, data = regression_data)   ### change phylop score category here
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(cassette_sums_heatmap)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = uniqueness_values, y = Cassette_Exon_Score)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(regression_data$uniqueness_values, na.rm = TRUE) * 0.6, y = max(regression_data$Cassette_Exon_Score, na.rm = TRUE) * 0.9, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Uniqueness Values (", i, ") for Cassette Exon"), y = "Mean PhyloP Score of Cassette Exon", title = paste0("Linear Regression: ", i, " Uniqueness Values vs PhyloP Score")) +
  theme_minimal()

print(p)
ggsave(filename = paste0("uniqueness_regression_plots/", i, "_regression_plot.png"), plot = p, width = 10, height = 6)
}

# Create a zip file of the plots:
zip(zipfile = "uniqueness_regression_plots.zip", files = "uniqueness_regression_plots")

# Optional: Remove the directory after zipping
unlink("uniqueness_regression_plots", recursive = TRUE)


# Now, instead of plotting uniqueness values, plot the individual % usage values of each cassette exon and correlate them with phylop scores:

regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)
temp_cassette <- mega_cassette[match(cassette_phylop$AS_event_ID, mega_cassette$AS_event_ID),]

# Create a directory to store plots and initialize empty list to store them in R:
dir.create("usage_regression_plots", showWarnings = FALSE)
plot_list <- list()

for (i in singletypes) {
  usage_values <- temp_cassette[[i]]
  usage_values <- as.numeric(gsub("%", "", usage_values))
  regression_data <- data.frame(usage_values, PhyloP_score = cassette_phylop$Cassette_Exon_Score)    ### change phylop score category here
  regression_data <- regression_data[complete.cases(regression_data),]
  regression_data <- regression_data[regression_data$PhyloP_score != "No Data",]
  regression_data$PhyloP_score <- as.numeric(regression_data$PhyloP_score)
  
  lm_model <- lm(PhyloP_score ~ usage_values, data = regression_data)
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(temp_cassette)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = as.numeric(usage_values), y = PhyloP_score)) +
  geom_point(color = "darkblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE, na.rm = TRUE) +
  annotate("text", x = 60, y = 5, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Percent Usage Values (", i, ") for Cassette Exon"), 
       y = "Mean PhyloP Score of Cassette Exon", 
       title = paste0("Linear Regression: ", i, " Percent Usage Values vs PhyloP Score")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_minimal()

print(p)
ggsave(filename = paste0("usage_regression_plots/", i, "_regression_plot.png"), plot = p, width = 10, height = 6)
}

# Create a zip file of the plots:
zip(zipfile = "usage_regression_plots.zip", files = "usage_regression_plots")

# Optional: Remove the directory after zipping
unlink("usage_regression_plots", recursive = TRUE)
