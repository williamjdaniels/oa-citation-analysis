# ==============================================================================
# Script 2: Analysis and Visualization
# Citation Analysis of Open Access Articles
# 
# Input:  min_data.csv, scaling_params.csv (from Script 1)
# Output: Figures, tables, model objects, and manuscript-ready exports
# ==============================================================================

# Load required packages
library(xtable)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(stargazer)
library(MASS)
library(emmeans)
library(lmtest)
library(DHARMa)
library(glmmTMB)
library(splines)
library(MuMIn)
library(performance)
library(scales)
library(patchwork)
library(car)
library(WeightIt)
library(cobalt)
library(broom)
library(dplyr)
library(knitr)
library(flextable)

# Set theme for all plots
theme_set(theme_classic() + 
            theme(text = element_text(size = 12),
                  axis.text = element_text(size = 10),
                  legend.text = element_text(size = 10)))

# Define color palette for OA types
oa_colors <- c("Bronze" = "#CD7F32", 
               "Closed Access" = "#000000", 
               "Gold" = "#DAA520", 
               "Green" = "#228B22")

# ==============================================================================
# Load processed data
# ==============================================================================
min_data <- data.table::fread("min_data.csv")

# Restore factor levels
min_data$`Open Access Stratified` <- factor(min_data$`Open Access Stratified`, 
                                            levels = c("Bronze", "Closed Access", "Green", "Gold"))
min_data$`Journal Title` <- factor(min_data$`Journal Title`)

# Load scaling parameters
scaling_params <- read.csv("scaling_params.csv")
author_mean <- scaling_params$mean[scaling_params$variable == "author_count"]
author_sd   <- scaling_params$sd[scaling_params$variable == "author_count"]
AIS_mean    <- scaling_params$mean[scaling_params$variable == "AIS"]
AIS_sd      <- scaling_params$sd[scaling_params$variable == "AIS"]
year_mean   <- scaling_params$mean[scaling_params$variable == "Publication Year"]
year_sd     <- scaling_params$sd[scaling_params$variable == "Publication Year"]

# ==============================================================================
# Descriptive Statistics
# ==============================================================================
descriptive_stats <- min_data %>%
  group_by(`Open Access Stratified`) %>%
  summarise(
    n = n(),
    mean_citations = mean(`Times Cited, All Databases`),
    median_citations = median(`Times Cited, All Databases`),
    sd_citations = sd(`Times Cited, All Databases`),
    iqr_citations = IQR(`Times Cited, All Databases`),
    zero_citations = sum(`Times Cited, All Databases` == 0),
    percent_zero = round(zero_citations / n * 100, 1)
  )

print("Descriptive Statistics by Open Access Type:")
print(descriptive_stats)

# ==============================================================================
# Generate Weights (IPTW)
# ==============================================================================
W.out <- weightit(`Open Access Stratified` ~ 
                    pub_year_scaled + 
                    auth_count_scaled + 
                    AIS_scaled + 
                    `Document Type`,
                  data = min_data, 
                  method = "ebal",
                  estimand = "ATE") 

# Check Balance
bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .1))
summary(W.out)

# Generate a Love Plot for the manuscript
Love <- love.plot(W.out, 
          thresholds = c(m = .1), 
          var.order = "unadjusted",
          abs = TRUE,
          line = TRUE, 
          stars = "raw") + 
  theme_classic() +
  labs(title = "Covariate Balance before and after weighting")
ggsave("Love_plot.tiff", plot = Love, width = 8, height = 8, dpi = 300)

# ==============================================================================
# Fit the final model
# ==============================================================================
min_data$weights <- W.out$weights

# final_weighted_model <- glmer.nb(`Times Cited, All Databases` ~ 
#                                    `Open Access Stratified` +
#                                    auth_count_scaled +
#                                    AIS_scaled +
#                                    pub_year_scaled +
#                                    `Document Type` +
#                                    pub_year_scaled:`Open Access Stratified` +
#                                    (pub_year_scaled | `Journal Title`),
#                                  data = min_data,
#                                  weights = weights, 
#                                  control = glmerControl(optimizer = "bobyqa",
#                                                         optCtrl = list(maxfun = 2e8)))
# final_weighted_model <- glmmTMB::glmmTMB(`Times Cited, All Databases` ~ 
#                                   `Open Access Stratified` +
#                                   auth_count_scaled +
#                                   AIS_scaled +
#                                   pub_year_scaled +
#                                   `Document Type` +
#                                   pub_year_scaled:`Open Access Stratified` +
#                                   (pub_year_scaled | `Journal Title`), # Random slope kept!
#                                 data = min_data,
#                                 weights = weights, 
#                                 family = nbinom2)
# model <- final_weighted_model
# saveRDS(final_weighted_model, "final_weighted_model.rds")

model <- readRDS("final_weighted_model.rds")
# Print summary to check results
summary(model)

# # check poisson
model_poisson <- glmmTMB::glmmTMB(`Times Cited, All Databases` ~ 
                                           `Open Access Stratified` +
                                           auth_count_scaled +
                                           AIS_scaled +
                                           pub_year_scaled +
                                           `Document Type` +
                                           pub_year_scaled:`Open Access Stratified` +
                                           (pub_year_scaled | `Journal Title`), # Random slope kept!
                                         data = min_data,
                                         weights = weights, 
                                         family = poisson)
# 
# # Then compare it to your negative binomial model
piosbinom <- anova(model_poisson, model)
piosbinom
summary(piosbinom)

# ==============================================================================
# Model diagnostics
# ==============================================================================
cat("\n=== Model Diagnostics ===\n")
overdisp_result <- check_overdispersion(model)
print(overdisp_result)

# Variance explained
r2_values <- r.squaredGLMM(model)
cat("\nVariance Explained (R²):\n")
cat(sprintf("Marginal R² (fixed effects): %.3f\n", r2_values[1]))
cat(sprintf("Conditional R² (fixed + random): %.3f\n", r2_values[2]))

# ICC for journal effects
# icc_value <- icc(model)
# cat(sprintf("\nIntraclass Correlation (Journal): %.3f\n", icc_value$ICC_adjusted))

# Extract model summary
model_summary <- summary(model)
print(model_summary)

# ==============================================================================
# Figure 1: Boxplot of citations by OA type
# ==============================================================================
fig1 <- ggplot(min_data, aes(x = `Open Access Stratified`, 
                             y = `Times Cited, All Databases`,
                             fill = `Open Access Stratified`)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1000")) +
  scale_fill_manual(values = oa_colors, guide = "none") +
  labs(x = "Open Access Type",
       y = "Citations (log scale)",
       title = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = list(c("Bronze", "Closed Access"),
                                        c("Bronze", "Gold"),
                                        c("Bronze", "Green"),
                                        c("Closed Access", "Gold"),
                                        c("Closed Access", "Green"),
                                        c("Gold", "Green")),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE)

# Add sample sizes to x-axis labels
n_labels <- descriptive_stats %>%
  mutate(label = paste0(`Open Access Stratified`, "\n(n=", scales::comma(n), ")"))

fig1 <- fig1 + scale_x_discrete(labels = n_labels$label)
ggsave("Access Type by Citations.tiff", fig1, width = 10, height = 10, dpi = 300)

# ==============================================================================
# Figure 2: Citations over time by OA type
# ==============================================================================
years <- sort(unique(min_data$`Publication Year`))
access_types <- levels(min_data$`Open Access Stratified`)
document_type <- (min_data$`Document Type`)

pred_grid_year <- expand.grid(
  `Publication Year` = years,
  `Open Access Stratified` = access_types,
  auth_count_scaled = 0,
  AIS_scaled = 0,
  `Document Type` = document_type,
  weights = 0
)

pred_grid_year$pub_year_scaled <- (pred_grid_year$`Publication Year` - year_mean) / year_sd
pred_grid_year$pub_year_scaled_sq <- pred_grid_year$pub_year_scaled^2

# Get predictions with confidence intervals
preds_year <- predict(model, newdata = pred_grid_year, type = "response", re.form = NA)
# preds_se_year <- predict(model, newdata = pred_grid_year, type = "link", re.form = NA, se.fit = TRUE)

# 1. Get the covariance matrix of the fixed effects
V <- vcov(model)$cond

# 2. Create the design matrix for your prediction grid
X <- model.matrix(~ `Open Access Stratified` + 
                    auth_count_scaled + 
                    AIS_scaled + 
                    pub_year_scaled + 
                    `Document Type` + 
                    pub_year_scaled:`Open Access Stratified`, 
                  data = pred_grid_year)

# 3. Calculate Standard Errors manually
se_link <- sqrt(rowSums((X %*% V) * X))

# 4. Get the predictions on the link scale (log scale)
pred_link <- X %*% fixef(model)$cond

# 5. Combine into your plotting object
plot_data_year <- pred_grid_year
plot_data_year$Citations <- exp(pred_link)
plot_data_year$Lower     <- exp(pred_link - 1.96 * se_link)
plot_data_year$Upper     <- exp(pred_link + 1.96 * se_link)

# plot_data_year <- pred_grid_year
# plot_data_year$Citations <- preds_year
# plot_data_year$Lower <- exp(preds_se_year$fit - 1.96 * preds_se_year$se.fit)
# plot_data_year$Upper <- exp(preds_se_year$fit + 1.96 * preds_se_year$se.fit)

YearPlot <- ggplot(plot_data_year, 
                   aes(color = `Open Access Stratified`,
                       fill = `Open Access Stratified`,
                       x = `Publication Year`,
                       y = Citations,
                       ymin = Lower, 
                       ymax = Upper)) +
  geom_pointrange(aes(shape = `Open Access Stratified`), 
                  position = position_dodge(width = 0.3), 
                  show.legend = TRUE) +
  geom_line(show.legend = FALSE) +
  labs(y = "Citations", x = "Year",
       title = "Relationship between year, access type, and citations") +
  theme_classic() +
  facet_wrap(~ `Document Type`, scales = "free_y") +
  scale_color_manual(values = c("#b08d57", "#000000", "#009E73", "#E69F00"), 
                     name = "Open Access Stratified") +
  scale_x_continuous(breaks = seq(min(plot_data_year$`Publication Year`), 
                                  max(plot_data_year$`Publication Year`), 
                                  by = 1)) +
  scale_shape_manual(values = c(15:18), name = "Open Access Stratified") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.9, 0.865))
ggsave("YearPlot.png", YearPlot, width = 12, height = 6, dpi = 300)

# ==============================================================================
# Figure 3: Citations by author count and OA type
# ==============================================================================
author_counts <- c(1, 2, 3, 5, 7, 10, 15)

pred_grid_author <- expand.grid(
  author_count = author_counts,
  `Open Access Stratified` = access_types,
  pub_year_scaled = 0,
  AIS_scaled = 0,
  `Document Type` = document_type,
  weights = 0
)

pred_grid_author$auth_count_scaled <- (pred_grid_author$author_count - author_mean) / author_sd
pred_grid_author$pub_year_scaled_sq <- 0

# Get predictions
preds_author <- predict(model, newdata = pred_grid_author, type = "response", re.form = NA)
# preds_se_author <- predict(model, newdata = pred_grid_author, 
#                            type = "link", 
#                            re.form = NA, 
#                            se.fit = TRUE)

# 2. Create the design matrix for the AUTHOR prediction grid
X_auth <- model.matrix(~ `Open Access Stratified` + 
                         auth_count_scaled + 
                         AIS_scaled + 
                         pub_year_scaled + 
                         `Document Type` + 
                         pub_year_scaled:`Open Access Stratified`, 
                       data = pred_grid_author)
se_link_author <- sqrt(rowSums((X_auth %*% V) * X_auth))

# 4. Get the predictions on the link scale (log scale)
pred_link_author <- X_auth %*% fixef(model)$cond

# 5. Combine into your plotting object
plot_data_author <- pred_grid_author
plot_data_author$Citations <- exp(pred_link_author)
plot_data_author$Lower     <- exp(pred_link_author - 1.96 * se_link_author)
plot_data_author$Upper     <- exp(pred_link_author + 1.96 * se_link_author)

AuthorPlot <- ggplot(plot_data_author, 
                     aes(color = `Open Access Stratified`,
                         fill = `Open Access Stratified`,
                         x = author_count,
                         y = Citations,
                         ymin = Lower, 
                         ymax = Upper)) +
  geom_pointrange(aes(shape = `Open Access Stratified`), 
                  position = position_dodge(width = 0.3), 
                  show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ `Document Type`, scales = "free_y") +
  labs(y = "Citations", x = "Number of Authors",
       title = "Relationship between author count, access type, and citations") +
  theme_classic() +
  scale_color_manual(values = c("#b08d57", "#000000", "#009E73", "#E69F00"), 
                     name = "Open Access Stratified") +
  scale_shape_manual(values = c(15:18), name = "Open Access Stratified") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.8, 0.8))

# Print the plot
print(AuthorPlot)
ggsave("AuthorPlot.tiff", AuthorPlot, width = 10, height = 6, dpi = 300)

# ==============================================================================
# Supplementary Figure: Pairwise comparisons
# ==============================================================================
emm <- emmeans(model, ~ `Open Access Stratified`, type = "response")
pairs_result <- pairs(emm, adjust = "tukey")

# Convert to data frame and calculate confidence intervals
pairs_df <- as.data.frame(pairs_result)

# Check column names
print("Column names in pairs_df:")
print(names(pairs_df))

# Calculate confidence intervals based on available columns
if("SE" %in% names(pairs_df)) {
  pairs_df$lower.CL <- pairs_df$ratio * exp(-1.96 * pairs_df$SE / pairs_df$ratio)
  pairs_df$upper.CL <- pairs_df$ratio * exp(1.96 * pairs_df$SE / pairs_df$ratio)
} else if("lower.CL" %in% names(pairs_df)) {
  pairs_df$lower.CL <- pairs_df$lower.CL
  pairs_df$upper.CL <- pairs_df$upper.CL
} else {
  pairs_ci <- confint(pairs_result)
  pairs_df <- cbind(pairs_df, pairs_ci)
}

pairs_df$significant <- ifelse(pairs_df$p.value < 0.05, "p < 0.05", "p ≥ 0.05")

supp_fig <- ggplot(pairs_df, aes(x = ratio, y = contrast, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("p < 0.05" = "red", "p ≥ 0.05" = "black"),
                     name = "Significance") +
  labs(x = "Citation Rate Ratio (95% CI)",
       y = "Pairwise Comparison",
       title = "Pairwise Comparisons of Citation Rates Between Open Access Types")
ggsave("citation_analysis_supplementary.tiff", supp_fig, width = 10, height = 6, dpi = 300)

# ==============================================================================
# Key Statistics for Manuscript
# ==============================================================================
cat("\n=== Key Statistics for Manuscript ===\n")
cat(sprintf("Total articles analyzed: %d\n", nrow(min_data)))
cat(sprintf("Number of journals: %d\n", n_distinct(min_data$`Journal Title`)))
cat(sprintf("Year range: %d-%d\n", min(min_data$`Publication Year`), max(min_data$`Publication Year`)))
cat(sprintf("Mean authors per article: %.1f (SD = %.1f)\n", author_mean, author_sd))
cat(sprintf("Percentage with zero citations: %.1f%%\n", 
            mean(min_data$`Times Cited, All Databases` == 0) * 100))

## remove authors for PLOS publication
min_data <- dplyr::select(min_data, !Authors)
# Save processed data for reproducibility
write.csv(min_data, "processed_citation_data.csv", row.names = FALSE)

# Save model object
saveRDS(model, "citation_model.rds")

# ==============================================================================
# Zero-Inflated Model Comparison
# ==============================================================================
zinb_model <- glmmTMB::glmmTMB(
  `Times Cited, All Databases` ~ 
    `Open Access Stratified` +
    auth_count_scaled +
    AIS_scaled +
    pub_year_scaled +
    `Document Type` +
    pub_year_scaled:`Open Access Stratified` +
    (pub_year_scaled | `Journal Title`), 
  data = min_data,
  weights = weights, 
  family = nbinom2,
  ziformula = ~1
)

# Compare the models
anova(model, zinb_model)
print("Analysis complete. Figures and data saved.")

# ==============================================================================
# ANOVA Table
# ==============================================================================
# Get Type II Anova table (recommended for models with interactions)
anova_table <- Anova(model, type = "II")
print(anova_table)

# Alternative: Use Type III if you prefer (requires setting contrasts)
# options(contrasts = c("contr.sum", "contr.poly"))
# anova_table_III <- Anova(model, type = "III")

# Format the output as a nice table
tidy_anova <- tidy(anova_table)

# Or use drop1 for a likelihood ratio test approach
drop1_table <- drop1(model, test = "Chisq")
print(drop1_table)

# For automatic extraction from your anova_table object:
anova_summary <- data.frame(
  Variable = rownames(anova_table),
  Chi_square = round(anova_table$Chisq, 2),
  df = anova_table$Df,
  p_value = ifelse(anova_table$`Pr(>Chisq)` < 0.001, 
                   "< 0.001", 
                   sprintf("%.3f", anova_table$`Pr(>Chisq)`))
)

# Clean up variable names for publication
anova_summary$Variable <- c("Open Access Type",
                            "Author Count",
                            "Article Influence Score", 
                            "Publication Year",
                            "Document Type",            
                            "Open Access Type × Year")
# Create publication-ready table
print(anova_summary)

# For LaTeX output
latex_table <- xtable(anova_summary,
                      caption = "Type II Wald Chi-square Tests for Fixed Effects",
                      label = "tab:anova",
                      align = c("l", "l", "r", "r", "r"))
names(latex_table) <- c("Variable", "$\\chi^2$", "df", "p-value")
print(latex_table, 
      include.rownames = FALSE,
      sanitize.text.function = function(x) x)

# For Word/HTML output using flextable
flex_table <- flextable(anova_summary)
flex_table <- set_header_labels(flex_table, 
                                Variable = "Variable",
                                Chi_square = "χ²",
                                df = "df",
                                p_value = "p-value")
flex_table <- theme_vanilla(flex_table)
flex_table <- autofit(flex_table)
save_as_docx(flex_table, path ="anova_table.docx")

# Save as CSV for easy import into Word/Excel
write.csv(anova_summary, "anova_results.csv", row.names = FALSE)

# ==============================================================================
# FINAL EXPORT SECTION: GENERATE MANUSCRIPT TABLES
# ==============================================================================

# 1. Export Table 1: Descriptive Statistics
write.csv(descriptive_stats, "Table1_Descriptive_Statistics.csv", row.names = FALSE)
cat("\n[Exported] Table 1 saved as 'Table1_Descriptive_Statistics.csv'")

# 2. Export Supplementary Table: Covariate Balance (Weighting Diagnostics)
balance_object <- bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .1))
balance_table <- balance_object$Balance.Across.Pairs

# Make row names (variable names) a proper column
balance_table$Variable <- row.names(balance_table)

# Reorder columns to put Variable first
balance_table <- balance_table %>% dplyr::select(Variable, everything())

write.csv(balance_table, "Supplement_Weighting_Balance.csv", row.names = FALSE)
cat("\n[Exported] Weighting balance saved as 'Supplement_Weighting_Balance.csv'")

# 3. Export Model Diagnostics (R-squared & Overdispersion)
diag_df <- data.frame(
  Metric = c("Marginal R2 (Fixed Effects)", 
             "Conditional R2 (Fixed + Random)", 
             "Overdispersion Ratio", 
             "Overdispersion p-value"),
  Value = c(r2_values[1], 
            r2_values[2], 
            overdisp_result$dispersion_ratio, 
            overdisp_result$p_value)
)

write.csv(diag_df, "Model_Diagnostics_Summary.csv", row.names = FALSE)
cat("\n[Exported] Model diagnostics saved as 'Model_Diagnostics_Summary.csv'")

# 4. Export Table 2: Full Model Estimates (Incidence Rate Ratios)
model_estimates <- broom.mixed::tidy(model, conf.int = TRUE, exponentiate = TRUE, effects = "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(
    term = gsub("Open Access Stratified", "OA: ", term),
    term = gsub("Document Type", "Type: ", term),
    across(c(estimate, conf.low, conf.high), ~round(., 3)),
    p.value = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  rename(
    Variable = term,
    IRR = estimate,
    CI_Lower = conf.low,
    CI_Upper = conf.high,
    P_Value = p.value
  )

write.csv(model_estimates, "Table2_Model_Estimates_IRR.csv", row.names = FALSE)
cat("\n[Exported] Model coefficients (IRRs) saved as 'Table2_Model_Estimates_IRR.csv'")

# 5. Export Table 3: Pairwise Comparisons (Post-Hoc Tests)
clean_pairs <- pairs_df %>%
  dplyr::select(contrast, ratio, p.value, lower.CL, upper.CL) %>%
  mutate(
    ratio = round(ratio, 3),
    lower.CL = round(lower.CL, 3),
    upper.CL = round(upper.CL, 3),
    p.value = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  rename(
    Comparison = contrast,
    Rate_Ratio = ratio,
    CI_Lower = lower.CL,
    CI_Upper = upper.CL,
    Adj_P_Value = p.value
  )

write.csv(clean_pairs, "Table3_Pairwise_Comparisons.csv", row.names = FALSE)
cat("\n[Exported] Pairwise comparisons saved as 'Table3_Pairwise_Comparisons.csv'")

# This gives you the table with Estimate (Beta), Std. Error, z value, and Pr(>|z|)
summary(model)$coefficients$cond

# Or for a cleaner version using the 'broom' package you already loaded:
broom.mixed::tidy(model, effects = "fixed")

emm_year <- emmeans(model, ~ `Open Access Stratified` | pub_year_scaled,
                     at = list(pub_year_scaled = (2014:2019 - year_mean) / year_sd),
                     type = "response")
pairs(emm_year, adjust = "tukey")

## 
cat("\n\n=== All Exports Complete ===")