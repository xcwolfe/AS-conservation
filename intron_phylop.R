library(tidyverse)

# Once you have run the Unix script average_phylop_score_for_IR_events.txt on HPCC, import your average scores for each IR event here:

intron_phylop <- read_tsv("phyloP_scores_combined.txt")
intron_phylop <- as.data.frame(intron_phylop)
intron_phylop$Average_Score <- as.numeric(intron_phylop$Average_Score)   # convert the Average_Score column to numeric
col_name <- "PhyloP score"
colnames(intron_phylop)[5] <- col_name

# Correlation of uniqueness index values with PhyloP score:

intron_sums_stats <- intron_sums_heatmap %>%
  rowwise() %>%
  summarize(AS_event_ID = AS_event_ID, Min_Value = min(c_across(-c(AS_event_ID, shift_or_preserve)), na.rm = TRUE), Max_Value = max(c_across(-c(AS_event_ID, shift_or_preserve)), na.rm = TRUE), shift_or_preserve = first(shift_or_preserve))

# "bin" min and max values into 5 separate percentiles:
intron_sums_stats <- intron_sums_stats %>%
  mutate(across(c(Min_Value, Max_Value), ~ cut(.x, quantile(.x, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                                               labels = c("0-20", "21-40", "41-60", "61-80", "81-100"), include.lowest = TRUE), 
                .names = "{.col}_Percentile"))

combined_data <- merge(intron_phylop, intron_sums_stats, by = "AS_event_ID")
combined_data <- combined_data[,!names(combined_data) %in% c("shift_or_preserve")]
combined_data <- combined_data[!is.na(combined_data$`PhyloP score`),]

lm_model <- lm((`Max_Value`) ~ `PhyloP score`, data = combined_data)
model_summary <- summary(lm_model)

slope <- coef(lm_model)[2]
intercept <- coef(lm_model)[1]
r_squared <- model_summary$adj.r.squared
p_value <- coef(model_summary)[2,4]

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nAdj. R² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

ggplot(combined_data, aes(x = Max_Value, y = `PhyloP score`)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(combined_data$Max_Value, na.rm = TRUE) * 0.7, y = max(combined_data$`PhyloP score`, na.rm = TRUE) * 0.9, label = custom_label, size = 5, color = "black", hjust = 0) +
  labs(x = "Max Uniqueness Value", y = "Mean PhyloP Score", title = "Linear Regression: Max Uniqueness Value vs Mean PhyloP Score for Each IR Event") +
  theme_minimal()

lm_model <- lm((`Min_Value`) ~ `PhyloP score`, data = combined_data)
model_summary <- summary(lm_model)

slope <- coef(lm_model)[2]
intercept <- coef(lm_model)[1]
r_squared <- model_summary$adj.r.squared
p_value <- coef(model_summary)[2,4]

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nAdj. R² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

ggplot(combined_data, aes(x = Min_Value, y = `PhyloP score`)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = min(combined_data$Min_Value, na.rm = TRUE) * 0.7, y = max(combined_data$`PhyloP score`, na.rm = TRUE) * 0.9, label = custom_label, size = 5, color = "black", hjust = 0) +
  labs(x = "Min Uniqueness Value", y = "Mean PhyloP Score", title = "Linear Regression: Min Uniqueness Value vs Mean PhyloP Score for Each IR Event")+
  theme_minimal()


# Do the same linear regression, but alter it for unique neuron types:

regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)

for (i in singletypes) {
  uniqueness_values <- intron_sums_heatmap[[i]]
  regression_data <- data.frame(uniqueness_values, PhyloP_score = intron_phylop$`PhyloP score`)
  regression_data <- regression_data[complete.cases(regression_data),]
  # eliminate zero-uniqueness values from graph/regression (optional, but may be helpful):
  regression_data <- regression_data[regression_data$uniqueness_values != 0,]
  
  lm_model <- lm(PhyloP_score ~ uniqueness_values, data = regression_data)
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(intron_sums_heatmap)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = uniqueness_values, y = PhyloP_score)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(regression_data$uniqueness_values, na.rm = TRUE) * 0.6, y = max(regression_data$PhyloP_score, na.rm = TRUE) * 0.9, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Uniqueness Values (", i, ") for Retained Intron"), y = "Mean PhyloP Score of Retained Intron", title = paste0("Linear Regression: ", i, " Uniqueness Values vs PhyloP Score")) +
  theme_minimal()

print(p)
}


# Now do a logistic regression with high, low, and mid uniqueness values:

combined_data <- merge(intron_phylop, intron_sums_stats, by = "AS_event_ID") %>%
  select(-shift_or_preserve) %>%
  filter(!is.na(`PhyloP score`)) %>%
  mutate(Max_Bin = factor(Max_Value_Percentile, levels = c("0-20", "21-40", "41-60", "61-80", "81-100")),
         Min_Bin = factor(Min_Value_Percentile, levels = c("0-20", "21-40", "41-60", "61-80", "81-100")))

logit_max <- glm(Max_Bin ~ `PhyloP score`, data = combined_data, family = binomial)
summary_max <- summary(logit_max)
logit_min <- glm(Min_Bin ~ `PhyloP score`, data = combined_data, family = binomial)
summary_min <- summary(logit_min)

p_value_max <- coef(summary_max)[2,4]    # p-value for PhyloP score in logit_max
p_value_min <- coef(summary_min)[2,4]

plot_logit <- function(model, data, xvar, yvar, title) {
  ggplot(data, aes(x = .data[[yvar]], y = .data[[xvar]], color = .data[[yvar]])) +
    geom_jitter(width = 0.35, height = 0.1, size = 1.8) +
    scale_x_discrete(labels = levels(data[[yvar]])) +  # Set x-axis labels to bin labels
    scale_color_manual(values = c("0-20" = "red", "21-40" = "red2", "41-60" = "red3", "61-80" = "maroon", "81-100" = "red4")) +
    labs(x = "Percentile Bin (higher percentile = higher uniqueness value)", y = xvar, title = title, color = "Percentile") +
    theme_minimal()
}

plot_max <- plot_logit(logit_max, combined_data, "PhyloP score", "Max_Bin", "Logistic Regression of IR events: PhyloP Score vs Max Uniqueness Value")
plot_min <- plot_logit(logit_min, combined_data, "PhyloP score", "Min_Bin", "Logistic Regression of IR events: PhyloP Score vs Min Uniqueness Value")

plot_max
plot_min


# Now, instead of plotting uniqueness values, plot the individual % usage values of each IR event and correlate them with phylop scores:

regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)
temp_intron <- mega_intron[match(intron_phylop$AS_event_ID, mega_intron$AS_event_ID),]

for (i in singletypes) {
  usage_values <- temp_intron[[i]]
  usage_values <- as.numeric(gsub("%", "", usage_values))
  regression_data <- data.frame(usage_values, PhyloP_score = intron_phylop$`PhyloP score`)
  regression_data <- regression_data[complete.cases(regression_data),]
  
  lm_model <- lm(PhyloP_score ~ usage_values, data = regression_data)
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(temp_intron)[i],
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
  labs(x = paste0("Percent Usage Values (", i, ") for Retained Intron"), 
       y = "Mean PhyloP Score of Retained Intron", 
       title = paste0("Linear Regression: ", i, " Percent Usage Values vs PhyloP Score")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_minimal()

print(p)
}


# Sort intron_phylop by highest to lowest phylop score:

intron_phylop <- intron_phylop[order(-intron_phylop$`PhyloP score`),]
intron_phylop$Gene <- mega_intron$Gene[match(intron_phylop$AS_event_ID, mega_intron$AS_event_ID)]
intron_phylop <- intron_phylop[,c("AS_event_ID","Chr","Start","End","Gene","PhyloP score")]
write.csv(intron_phylop, file = "intron_retention_phylop_scores.csv", row.names = F)

