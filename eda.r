# =============================================================================
# Chernobyl Radiation Data — Exploratory Data Analysis
# =============================================================================
# Author  : Ashitom Budhathoki
# Data    : chernobyl.csv
# Columns : PAYS (Country), Ville (City), Date, X, Y,
#           I 131 (Bq/m3), Cs 134 (Bq/m3), Cs 137 (Bq/m3)
# =============================================================================


# ── 0. Dependencies ───────────────────────────────────────────────────────────

required_packages <- c("tidyverse", "ggcorrplot", "scales", "patchwork", "ggridges")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

invisible(lapply(required_packages, install_if_missing))

library(tidyverse)
library(ggcorrplot)
library(scales)
library(patchwork)
library(ggridges)


# ── 1. Configuration ──────────────────────────────────────────────────────────

DATA_PATH   <- "chernobyl.csv" 
OUTPUT_DIR  <- "eda_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

ISOTOPES <- c("I_131", "Cs_134", "Cs_137")   # clean column names (set below)
PALETTE  <- c("#2196F3", "#FF5722", "#4CAF50") # blue, orange, green

theme_chernobyl <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 14, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 11, colour = "grey40"),
      axis.title    = element_text(size = 11),
      legend.position = "bottom",
      strip.text    = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(plot, filename, w = 10, h = 6) {
  path <- file.path(OUTPUT_DIR, filename)
  ggsave(path, plot = plot, width = w, height = h, dpi = 150)
  message("  Saved: ", path)
}


# ── 2. Data Ingestion & Cleaning ──────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("  CHERNOBYL EDA — Loading & Cleaning Data\n")
cat(strrep("=", 60), "\n\n")

raw <- read_csv(DATA_PATH, show_col_types = FALSE)

df <- raw |>
  rename(
    Country = PAYS,
    City    = Ville,
    I_131   = `I 131 (Bq/m3)`,
    Cs_134  = `Cs 134 (Bq/m3)`,
    Cs_137  = `Cs 137 (Bq/m3)`,
    Lon     = X,
    Lat     = Y
  ) |>
  mutate(
    across(c(I_131, Cs_134, Cs_137), ~ suppressWarnings(as.numeric(.))),
    Date    = as.Date(Date),
    Country = as.factor(Country),
    City    = as.factor(City)
  )

cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
cat("Columns          :", paste(names(df), collapse = ", "), "\n\n")


# ── 3. Summary Statistics ─────────────────────────────────────────────────────

cat(strrep("-", 60), "\n")
cat("  SECTION 1 — Summary Statistics\n")
cat(strrep("-", 60), "\n\n")

# 3a. Missing values
cat("Missing values:\n")
df |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "Column", values_to = "Missing") |>
  mutate(Pct = scales::percent(Missing / nrow(df))) |>
  filter(Missing > 0) |>
  print()

# 3b. Descriptive stats for isotopes
cat("\nDescriptive statistics (Bq/m³):\n")
df |>
  select(all_of(ISOTOPES)) |>
  pivot_longer(everything(), names_to = "Isotope", values_to = "Bq_m3") |>
  group_by(Isotope) |>
  summarise(
    n      = sum(!is.na(Bq_m3)),
    mean   = mean(Bq_m3, na.rm = TRUE),
    median = median(Bq_m3, na.rm = TRUE),
    sd     = sd(Bq_m3, na.rm = TRUE),
    min    = min(Bq_m3, na.rm = TRUE),
    max    = max(Bq_m3, na.rm = TRUE),
    q25    = quantile(Bq_m3, 0.25, na.rm = TRUE),
    q75    = quantile(Bq_m3, 0.75, na.rm = TRUE),
    iqr    = IQR(Bq_m3, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print(n = Inf)

# 3c. Country-wise mean
cat("\nCountry-wise mean radiation (Bq/m³):\n")
df |>
  group_by(Country) |>
  summarise(across(all_of(ISOTOPES), ~ mean(., na.rm = TRUE)), .groups = "drop") |>
  arrange(desc(I_131)) |>
  print(n = Inf)

# 3d. Top 10 cities by total mean radiation
cat("\nTop 10 cities by total mean radiation:\n")
df |>
  group_by(City) |>
  summarise(across(all_of(ISOTOPES), ~ mean(., na.rm = TRUE)), .groups = "drop") |>
  mutate(Total = I_131 + Cs_134 + Cs_137) |>
  slice_max(Total, n = 10) |>
  print()


# ── 4. Outlier Detection ──────────────────────────────────────────────────────

cat("\n", strrep("-", 60), "\n")
cat("  SECTION 2 — Outlier Detection\n")
cat(strrep("-", 60), "\n\n")

iqr_outliers <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  sum(x < (q[1] - 1.5 * iqr) | x > (q[2] + 1.5 * iqr), na.rm = TRUE)
}
zscore_outliers <- function(x, threshold = 3) {
  z <- abs((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  sum(z > threshold, na.rm = TRUE)
}

df |>
  select(all_of(ISOTOPES)) |>
  pivot_longer(everything(), names_to = "Isotope", values_to = "val") |>
  group_by(Isotope) |>
  summarise(
    IQR_outliers    = iqr_outliers(val),
    Zscore_outliers = zscore_outliers(val),
    .groups = "drop"
  ) |>
  print()


# ── 5. Correlation Analysis ───────────────────────────────────────────────────

cat("\n", strrep("-", 60), "\n")
cat("  SECTION 3 — Correlation Analysis\n")
cat(strrep("-", 60), "\n\n")

corr_mat <- df |>
  select(all_of(ISOTOPES)) |>
  cor(use = "pairwise.complete.obs")

cat("Pearson correlation matrix:\n")
print(round(corr_mat, 4))


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION A — Categorical / Sample Distribution Plots
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n", strrep("=", 60), "\n")
cat("  PLOTS — Section A: Sample Distributions\n")
cat(strrep("=", 60), "\n")

# A1. Samples by Country
p_a1 <- df |>
  count(Country, sort = TRUE) |>
  ggplot(aes(x = fct_reorder(Country, n), y = n)) +
  geom_col(fill = "#2196F3") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Sample Distribution by Country",
    subtitle = "Total number of radiation measurements per country",
    x = NULL, y = "Number of Samples"
  ) +
  theme_chernobyl()

save_plot(p_a1, "A1_country_distribution.png", w = 9, h = 5)

# A2. Top 15 Cities
p_a2 <- df |>
  count(City, sort = TRUE) |>
  slice_head(n = 15) |>
  ggplot(aes(x = fct_reorder(City, n), y = n)) +
  geom_col(fill = "#FF5722") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Sample Distribution by City (Top 15)",
    subtitle = "Cities with the most radiation measurement records",
    x = NULL, y = "Number of Samples"
  ) +
  theme_chernobyl()

save_plot(p_a2, "A2_city_distribution.png", w = 9, h = 5)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION B — Univariate Distribution Plots
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 60), "\n")
cat("  PLOTS — Section B: Univariate Distributions\n")
cat(strrep("=", 60), "\n")

df_long <- df |>
  select(all_of(ISOTOPES)) |>
  pivot_longer(everything(), names_to = "Isotope", values_to = "Bq_m3") |>
  filter(!is.na(Bq_m3))

# B1. Raw histogram (clipped)
p_b1 <- df_long |>
  filter(Bq_m3 <= 10) |>
  ggplot(aes(x = Bq_m3, fill = Isotope)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = PALETTE) +
  scale_x_continuous(labels = comma) +
  labs(
    title    = "Distribution of Radiation Measurements",
    subtitle = "Clipped to ≤ 10 Bq/m³ for readability",
    x = "Bq/m³", y = "Frequency", fill = "Isotope"
  ) +
  theme_chernobyl()

save_plot(p_b1, "B1_radiation_histogram.png")

# B2. Log-transformed histogram
p_b2 <- df_long |>
  ggplot(aes(x = log1p(Bq_m3), fill = Isotope)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = PALETTE) +
  labs(
    title    = "Log-Transformed Distribution of Radiation",
    subtitle = "log(1 + Bq/m³) — reveals full range of values",
    x = "log(1 + Bq/m³)", y = "Frequency", fill = "Isotope"
  ) +
  theme_chernobyl()

save_plot(p_b2, "B2_log_histogram.png")

# B3. Density (log scale)
p_b3 <- df_long |>
  filter(Bq_m3 > 0) |>
  ggplot(aes(x = log1p(Bq_m3), fill = Isotope, colour = Isotope)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values = PALETTE) +
  scale_colour_manual(values = PALETTE) +
  labs(
    title    = "Log-Transformed Density Plot",
    subtitle = "Smoothed kernel density of log(1 + Bq/m³)",
    x = "log(1 + Bq/m³)", y = "Density", fill = "Isotope", colour = "Isotope"
  ) +
  theme_chernobyl()

save_plot(p_b3, "B3_log_density.png")

# B4. Box plot
p_b4 <- df_long |>
  filter(Bq_m3 <= 5) |>
  ggplot(aes(x = Isotope, y = Bq_m3, fill = Isotope)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4, width = 0.5) +
  scale_fill_manual(values = PALETTE) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Box Plot of Radiation Measurements",
    subtitle = "Clipped to ≤ 5 Bq/m³; dots are outliers",
    x = NULL, y = "Bq/m³"
  ) +
  theme_chernobyl() +
  theme(legend.position = "none")

save_plot(p_b4, "B4_boxplot.png", w = 7, h = 5)

# B5. Q-Q Plot (I-131)
p_b5 <- df |>
  filter(!is.na(I_131), I_131 > 0) |>
  ggplot(aes(sample = log1p(I_131))) +
  stat_qq(colour = PALETTE[1], alpha = 0.5, size = 0.8) +
  stat_qq_line(colour = "black", linewidth = 0.7) +
  labs(
    title    = "Q-Q Plot: I-131 (log scale)",
    subtitle = "Assessing normality of log(1 + I-131)",
    x = "Theoretical Quantiles", y = "Sample Quantiles"
  ) +
  theme_chernobyl()

save_plot(p_b5, "B5_qq_plot_i131.png", w = 7, h = 5)

# B6. Ridge plot — isotope distributions by country
p_b6 <- df |>
  select(Country, all_of(ISOTOPES)) |>
  pivot_longer(all_of(ISOTOPES), names_to = "Isotope", values_to = "Bq_m3") |>
  filter(!is.na(Bq_m3), Bq_m3 > 0) |>
  ggplot(aes(x = log1p(Bq_m3), y = Country, fill = Isotope)) +
  geom_density_ridges(alpha = 0.55, rel_min_height = 0.01, scale = 0.9) +
  scale_fill_manual(values = PALETTE) +
  facet_wrap(~Isotope, ncol = 3) +
  labs(
    title    = "Ridge Plot: log(1 + Bq/m³) by Country",
    subtitle = "Distribution shape of each isotope across countries",
    x = "log(1 + Bq/m³)", y = NULL
  ) +
  theme_chernobyl() +
  theme(legend.position = "none")

save_plot(p_b6, "B6_ridge_by_country.png", w = 14, h = 7)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION C — Country & City Comparative Plots
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 60), "\n")
cat("  PLOTS — Section C: Country & City Comparisons\n")
cat(strrep("=", 60), "\n")

# C1. Mean radiation by country (grouped bar)
p_c1 <- df |>
  group_by(Country) |>
  summarise(across(all_of(ISOTOPES), ~ mean(., na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(ISOTOPES), names_to = "Isotope", values_to = "Mean_Bq") |>
  ggplot(aes(x = fct_reorder(Country, Mean_Bq, .fun = sum), y = Mean_Bq, fill = Isotope)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = PALETTE) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Mean Radiation by Country",
    subtitle = "Average Bq/m³ per isotope grouped by country",
    x = NULL, y = "Mean Bq/m³", fill = "Isotope"
  ) +
  theme_chernobyl()

save_plot(p_c1, "C1_mean_by_country.png", w = 10, h = 6)

# C2. Top 10 cities by mean radiation (horizontal stacked bar)
p_c2 <- df |>
  group_by(City) |>
  summarise(across(all_of(ISOTOPES), ~ mean(., na.rm = TRUE)), .groups = "drop") |>
  mutate(Total = I_131 + Cs_134 + Cs_137) |>
  slice_max(Total, n = 10) |>
  select(-Total) |>
  pivot_longer(all_of(ISOTOPES), names_to = "Isotope", values_to = "Mean_Bq") |>
  ggplot(aes(x = fct_reorder(City, Mean_Bq, .fun = sum), y = Mean_Bq, fill = Isotope)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual(values = PALETTE) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Top 10 Cities by Mean Radiation",
    subtitle = "Stacked average Bq/m³ across all three isotopes",
    x = NULL, y = "Mean Bq/m³", fill = "Isotope"
  ) +
  theme_chernobyl()

save_plot(p_c2, "C2_top10_cities.png", w = 10, h = 6)

# C3. Violin plot by country
p_c3 <- df_long |>
  filter(Bq_m3 > 0, Bq_m3 <= 5) |>
  ggplot(aes(x = Country, y = Bq_m3, fill = Country)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title    = "Violin Plot of Radiation by Country",
    subtitle = "Combined distribution of all three isotopes (≤ 5 Bq/m³)",
    x = NULL, y = "Bq/m³"
  ) +
  theme_chernobyl() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(p_c3, "C3_violin_by_country.png", w = 12, h = 6)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION D — Bivariate / Correlation Plots
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 60), "\n")
cat("  PLOTS — Section D: Bivariate & Correlation Analysis\n")
cat(strrep("=", 60), "\n")

# D1. Correlation heatmap
p_d1 <- ggcorrplot(
  corr_mat,
  method   = "square",
  type     = "upper",
  lab      = TRUE,
  lab_size = 4,
  colors   = c("#D32F2F", "white", "#1976D2"),
  title    = "Correlation Matrix (Pearson)",
  ggtheme  = theme_chernobyl()
) +
  labs(subtitle = "Upper triangle; values show Pearson r")

save_plot(p_d1, "D1_correlation_heatmap.png", w = 7, h = 6)

# D2. Scatter: I-131 vs Cs-137
p_d2 <- df |>
  filter(!is.na(I_131), !is.na(Cs_137)) |>
  ggplot(aes(x = I_131, y = Cs_137)) +
  geom_point(alpha = 0.4, colour = PALETTE[1], size = 1.2) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  scale_x_continuous(limits = c(0, quantile(df$I_131, 0.99, na.rm = TRUE))) +
  scale_y_continuous(limits = c(0, quantile(df$Cs_137, 0.99, na.rm = TRUE))) +
  labs(
    title    = "I-131 vs Cs-137",
    subtitle = "Clipped to 99th percentile; black line = OLS fit",
    x = "I-131 (Bq/m³)", y = "Cs-137 (Bq/m³)"
  ) +
  theme_chernobyl()

save_plot(p_d2, "D2_scatter_I131_vs_Cs137.png")

# D3. Scatter: Cs-134 vs Cs-137
p_d3 <- df |>
  filter(!is.na(Cs_134), !is.na(Cs_137)) |>
  ggplot(aes(x = Cs_134, y = Cs_137)) +
  geom_point(alpha = 0.4, colour = PALETTE[2], size = 1.2) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  scale_x_continuous(limits = c(0, quantile(df$Cs_134, 0.99, na.rm = TRUE))) +
  scale_y_continuous(limits = c(0, quantile(df$Cs_137, 0.99, na.rm = TRUE))) +
  labs(
    title    = "Cs-134 vs Cs-137",
    subtitle = "Clipped to 99th percentile; black line = OLS fit",
    x = "Cs-134 (Bq/m³)", y = "Cs-137 (Bq/m³)"
  ) +
  theme_chernobyl()

save_plot(p_d3, "D3_scatter_Cs134_vs_Cs137.png")

# D4. Pairs plot (log-scale)
pairs_data <- df |>
  select(all_of(ISOTOPES)) |>
  filter(complete.cases(pick(everything()))) |>
  mutate(across(everything(), log1p))

png(file.path(OUTPUT_DIR, "D4_pairs_plot.png"), width = 800, height = 800, res = 120)
pairs(
  pairs_data,
  labels    = ISOTOPES,
  col       = adjustcolor(PALETTE[1], alpha.f = 0.3),
  pch       = 19,
  cex       = 0.4,
  gap       = 0.5,
  main      = "Pairs Plot (log scale)"
)
dev.off()
message("  Saved: ", file.path(OUTPUT_DIR, "D4_pairs_plot.png"))


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION E — Geospatial Plots
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 60), "\n")
cat("  PLOTS — Section E: Geographic Distribution\n")
cat(strrep("=", 60), "\n")

df_geo <- df |>
  filter(!is.na(Lon), !is.na(Lat), !is.na(I_131), !is.na(Cs_137))

# E1. Geographic scatter — I-131
p_e1 <- df_geo |>
  ggplot(aes(x = Lon, y = Lat, colour = I_131)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_colour_gradient(low = "lightyellow", high = "red3", labels = comma,
                        name = "I-131 (Bq/m³)") +
  labs(
    title    = "Geographic Distribution — I-131",
    subtitle = "Point colour reflects radiation intensity",
    x = "Longitude", y = "Latitude"
  ) +
  theme_chernobyl() +
  theme(legend.position = "right")

save_plot(p_e1, "E1_geo_I131.png", w = 10, h = 7)

# E2. Geographic scatter — Cs-137
p_e2 <- df_geo |>
  ggplot(aes(x = Lon, y = Lat, colour = Cs_137)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_colour_gradient(low = "lightyellow", high = "darkblue", labels = comma,
                        name = "Cs-137 (Bq/m³)") +
  labs(
    title    = "Geographic Distribution — Cs-137",
    subtitle = "Point colour reflects radiation intensity",
    x = "Longitude", y = "Latitude"
  ) +
  theme_chernobyl() +
  theme(legend.position = "right")

save_plot(p_e2, "E2_geo_Cs137.png", w = 10, h = 7)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION F — Temporal Plot
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n", strrep("=", 60), "\n")
cat("  PLOTS — Section F: Temporal Overview\n")
cat(strrep("=", 60), "\n")

# F1. Weekly mean over time
p_f1 <- df |>
  filter(!is.na(Date)) |>
  mutate(Week = lubridate::floor_date(Date, "week")) |>
  group_by(Week) |>
  summarise(across(all_of(ISOTOPES), ~ mean(., na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(ISOTOPES), names_to = "Isotope", values_to = "Mean_Bq") |>
  ggplot(aes(x = Week, y = Mean_Bq, colour = Isotope)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = PALETTE) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Weekly Mean Radiation Over Time",
    subtitle = "Each point = weekly average Bq/m³ across all locations",
    x = NULL, y = "Mean Bq/m³", colour = "Isotope"
  ) +
  theme_chernobyl()

save_plot(p_f1, "F1_weekly_mean_over_time.png", w = 12, h = 5)


# ── 6. Final Summary ──────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("  EDA COMPLETE — Summary\n")
cat(strrep("=", 60), "\n\n")

cat(sprintf("  Total records    : %d\n",       nrow(df)))
cat(sprintf("  Countries        : %d\n",       n_distinct(df$Country)))
cat(sprintf("  Cities           : %d\n",       n_distinct(df$City)))
cat(sprintf("  Date range       : %s to %s\n", min(df$Date, na.rm = TRUE),
                                                max(df$Date, na.rm = TRUE)))
cat(sprintf("  Missing (isotopes): %d\n",      sum(is.na(df[ISOTOPES]))))
cat(sprintf("  Output directory : %s/\n\n",    OUTPUT_DIR))

all_plots <- list.files(OUTPUT_DIR, pattern = "\\.png$")
cat("  Plots saved:\n")
cat(paste0("    ", seq_along(all_plots), ". ", all_plots, "\n"), sep = "")

cat("\n", strrep("=", 60), "\n")
cat("  Done!\n")
cat(strrep("=", 60), "\n\n")
