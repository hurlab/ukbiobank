#!/usr/bin/env Rscript

################################################################################
# UK Biobank: Combine OLINK protein measurements with MNSI questionnaire responses
# MNSI scoring is performed inline. If the raw MNSI is NA due to missing items,
# keep the partial sum when the sum of available (non-NA) scored items is ≥ 3
################################################################################

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------

# Set working directory to this script's location (RStudio or command line)
set_wd_to_script <- function() {
  pth <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (rlang::is_true(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
      pth <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    }
  }
  if (is.null(pth) || !nzchar(pth)) {
    # Fallback for Rscript execution
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1L) pth <- sub("^--file=", "", file_arg)
  }
  if (!is.null(pth) && nzchar(pth)) setwd(dirname(normalizePath(pth)))
}
set_wd_to_script()

# General options
options(stringsAsFactors = FALSE)

# Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table, readr, dplyr, tidyr, stringr, tibble, purrr, glue, magrittr,
  ggplot2, VennDetail, rlang
)

# Output directories
out_dir <- file.path("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figs"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "data"),   showWarnings = FALSE, recursive = TRUE)

## -----------------------------------------------------------------------------
## Configuration: file paths
## -----------------------------------------------------------------------------
olink_path    <- file.path("rawdata", "plasma", "OLINK protein with name_instance0.csv")
question_path <- file.path("rawdata", "plasma", "Questionnaire_v2.csv")
# write combined CSV into the 'data' output directory
output_path   <- file.path(out_dir, "data", "olink_questionnaire_matched.csv")

## -----------------------------------------------------------------------------
## I/O checks
## -----------------------------------------------------------------------------
if (!file.exists(olink_path)) {
  stop("OLINK protein file not found: ", olink_path)
}
if (!file.exists(question_path)) {
  stop("Questionnaire file not found: ", question_path)
}

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Map questionnaire responses to 0/1, accepting "YES"/"NO" or numeric 0/1
map_response <- function(x) {
  if (is.numeric(x)) {
    return(ifelse(x %in% c(0, 1), as.numeric(x), NA_real()))
  }
  nx <- trimws(toupper(as.character(x)))
  out <- rep(NA_real_, length(nx))
  out[nx == "YES"] <- 1
  out[nx == "NO"]  <- 0
  as.numeric(out)
}

# Treat negative values as missing
clean_negatives <- function(v) {
  v[v < 0] <- NA
  v
}

## -----------------------------------------------------------------------------
## Load data
## -----------------------------------------------------------------------------

# OLINK
olink_raw <- read.csv(olink_path, check.names = FALSE)
colnames(olink_raw)[1] <- "eid"
protein_labels  <- colnames(olink_raw)
protein_symbols <- protein_labels
protein_symbols[-1] <- trimws(sub(";.*$", "", protein_labels[-1]))
colnames(olink_raw) <- make.unique(protein_symbols)
olink_raw$eid <- as.character(olink_raw$eid)

# Questionnaire v2
question_raw <- read.csv(question_path, check.names = FALSE)

# Standardize anchors
colnames(question_raw)[1] <- "record_id"

# MNSI code columns still start right after the first column and span 15 items
question_codes <- 120071:120085
question_cols  <- sprintf("f.%06d.0.0", question_codes)
stopifnot(length(question_codes) == 15)
colnames(question_raw)[2:(length(question_cols) + 1)] <- question_cols

# Two extra columns inserted right before ID in v2:
#   "Ever had diabetes (Type I or Type II)"
#   "Age at recruitment"  -> rename to 'age_at_recruitment'
if ("Age at recruitment" %in% names(question_raw)) {
  names(question_raw)[names(question_raw) == "Age at recruitment"] <- "age_at_recruitment"
}

# ID and sex handling (ID is no longer the last column; sex comes right after ID)
if (!"ID" %in% names(question_raw)) {
  stop("Expected 'ID' column not found in questionnaire file.")
}
question_raw$eid <- as.character(question_raw$ID)

# Normalize questionnaire responses to 0/1 for MNSI items
question_raw[question_cols] <- lapply(question_raw[question_cols], map_response)

# Sex column (1 = male, 0 = female). Keep as integer for downstream use.
if ("sex" %in% names(question_raw)) {
  question_raw$sex <- suppressWarnings(as.integer(question_raw$sex))
}

# Diabetes flag
# Prefer explicit v2 column if present; otherwise fall back to UKB 120007 logic or default 1
diab_col_v2 <- "Ever had diabetes (Type I or Type II)"
if (diab_col_v2 %in% names(question_raw)) {
  diab_mapped <- map_response(question_raw[[diab_col_v2]])
  question_raw$f.120007.0.0 <- diab_mapped
} else if (!"f.120007.0.0" %in% names(question_raw)) {
  question_raw$f.120007.0.0 <- 1L
}

# Drop duplicate questionnaire rows by eid, keep first
dup_idx <- duplicated(question_raw$eid)
if (any(dup_idx)) {
  warning("Duplicate questionnaire IDs detected; keeping first occurrence.")
  question_raw <- question_raw[!dup_idx, ]
}

## -----------------------------------------------------------------------------
## Match participants present in both datasets
## -----------------------------------------------------------------------------
common_ids <- sort(intersect(question_raw$eid, olink_raw$eid))
if (!length(common_ids)) stop("No participants found in both datasets.")

question_matched <- question_raw[match(common_ids, question_raw$eid), ]
olink_matched    <- olink_raw[match(common_ids, olink_raw$eid), ]

missing_rows <- is.na(question_matched$eid) | is.na(olink_matched$eid)
if (any(missing_rows)) {
  question_matched <- question_matched[!missing_rows, , drop = FALSE]
  olink_matched    <- olink_matched[!missing_rows, , drop = FALSE]
}

## -----------------------------------------------------------------------------
## MNSI scoring
## -----------------------------------------------------------------------------
# Items 120071 to 120085 are the questionnaire. Two are reversed: 120077, 120083.
# Two items are not scored in the standard self-admin MNSI: 120074 (Q4) and 120080 (Q10).
# Compute raw sum across scored items. For incomplete rows, keep partial sum if > 3.

mnsi_all_cols <- sprintf("f.%06d.0.0", 120071:120085)
present_cols  <- intersect(mnsi_all_cols, names(question_matched))
question_matched[present_cols] <- lapply(question_matched[present_cols], clean_negatives)

# Build reversed items
if ("f.120083.0.0" %in% names(question_matched)) {
  question_matched$f.120083.0.0_rev <- ifelse(
    question_matched$f.120083.0.0 == 1, 0,
    ifelse(question_matched$f.120083.0.0 == 0, 1, NA)
  )
} else {
  question_matched$f.120083.0.0_rev <- NA_real_
}

if ("f.120077.0.0" %in% names(question_matched)) {
  question_matched$f.120077.0.0_rev <- ifelse(
    question_matched$f.120077.0.0 == 1, 0,
    ifelse(question_matched$f.120077.0.0 == 0, 1, NA)
  )
} else {
  question_matched$f.120077.0.0_rev <- NA_real_
}

# Scored item set (exclude 120074 and 120080)
mnsi_items <- c(
  "f.120071.0.0",
  "f.120072.0.0",
  "f.120073.0.0",
  "f.120075.0.0",
  "f.120076.0.0",
  "f.120077.0.0_rev",
  "f.120078.0.0",
  "f.120079.0.0",
  "f.120081.0.0",
  "f.120082.0.0",
  "f.120083.0.0_rev",
  "f.120084.0.0",
  "f.120085.0.0"
)

## -----------------------------------------------------------------------------
## Initialize outputs and MNSI scoring (compute for all participants)
## -----------------------------------------------------------------------------

# Keep diabetes for reference, but do not gate MNSI on it
question_matched$diabetes <- ifelse(
  is.na(question_matched$f.120007.0.0), NA_integer_,
  ifelse(question_matched$f.120007.0.0 == 1, 1L, 0L)
)

question_matched$MNSI <- NA_real_
question_matched$pn   <- NA_integer_

# Build matrix of MNSI items for all participants
M <- question_matched[, mnsi_items, drop = FALSE]

# Sums and missingness
raw_sum           <- rowSums(M, na.rm = FALSE)  # NA if any NA present
partial_sum       <- rowSums(M, na.rm = TRUE)   # ignores NAs
non_missing_count <- rowSums(!is.na(M))
all_missing       <- non_missing_count == 0
has_missing       <- non_missing_count < ncol(M)

# Start with raw sums for complete cases
final_sum <- raw_sum

# For rows with any missing items:
# keep partial_sum if > 3, else set to NA
keep_partial <- has_missing & (partial_sum > 3)
drop_partial <- has_missing & !keep_partial

final_sum[keep_partial] <- partial_sum[keep_partial]
final_sum[drop_partial] <- NA_real_
final_sum[all_missing]  <- NA_real_

# Write back for all participants
question_matched$MNSI <- final_sum

# PN outcome: 1 if MNSI > 3, 0 if MNSI ≤ 3, NA if MNSI is NA
question_matched$pn <- ifelse(
  is.na(question_matched$MNSI), NA_integer_,
  ifelse(question_matched$MNSI > 3, 1L, 0L)
)

## -----------------------------------------------------------------------------
## Combine and write
## -----------------------------------------------------------------------------
combined <- cbind(question_matched, olink_matched[, -1, drop = FALSE])
write.csv(combined, output_path, row.names = FALSE)
message("Wrote matched protein and questionnaire data to ", output_path)


## -----------------------------------------------------------------------------
## Participant summary by neuropathy status
## -----------------------------------------------------------------------------

# Exclude records with missing diabetes, MNSI, or pn
valid_df <- combined %>%
  filter(!is.na(diabetes), !is.na(MNSI), !is.na(pn))

# Define groups
valid_df <- valid_df %>%
  mutate(group = case_when(
    diabetes == 1 & pn == 1 ~ "diabetes_DN",
    diabetes == 1 & pn == 0 ~ "diabetes_nonDN",
    TRUE ~ NA_character_
  ))

# Summarize group counts
group_counts <- valid_df %>%
  filter(!is.na(group)) %>%
  count(group, name = "n_participants")

# Print results
cat("\n================ Participant Counts ================\n")
print(group_counts)
cat("====================================================\n")

# Optionally save the summary table
write.csv(group_counts, file.path(out_dir, "tables", "group_counts.csv"), row.names = FALSE)
message("Saved participant group counts to output/tables/group_counts.csv")


## -----------------------------------------------------------------------------
## Differential protein analysis: proteins from column 28 onward (robust version)
## -----------------------------------------------------------------------------
library(broom)
library(purrr)
library(dplyr)
library(tidyr)

if (ncol(combined) < 28) stop("Expected protein markers to start at column 28, but ncol(combined) < 28.")

meta_keep   <- c("record_id", "eid", "diabetes", "MNSI", "pn", "sex", "age_at_recruitment")
protein_cols <- colnames(combined)[28:ncol(combined)]

# Build analysis dataframe and set group factor with stable reference
analysis_df <- combined %>%
  select(any_of(meta_keep), all_of(protein_cols)) %>%
  filter(!is.na(diabetes), !is.na(MNSI), !is.na(pn), !is.na(sex), !is.na(age_at_recruitment)) %>%
  mutate(group = case_when(
    diabetes == 1 & pn == 1 ~ "diabetes_DN",
    diabetes == 1 & pn == 0 ~ "diabetes_nonDN",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  mutate(
    group = factor(group, levels = c("diabetes_nonDN", "diabetes_DN")),
    sex   = as.integer(sex) # keep numeric 0/1
  )

# Long format
df_long <- analysis_df %>%
  pivot_longer(cols = all_of(protein_cols), names_to = "Protein", values_to = "NPX")

message("Running linear regression for ", length(unique(df_long$Protein)), " proteins...")

# Helper to fit one model and return tidy row for the group effect
fit_one <- function(dat) {
  mdl <- tryCatch(lm(NPX ~ group + age_at_recruitment + sex, data = dat), error = function(e) NULL)
  if (is.null(mdl)) {
    return(tibble(term = "groupdiabetes_DN", estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_))
  }
  tm <- tidy(mdl)
  row <- tm %>% filter(term == "groupdiabetes_DN") %>% slice_head(n = 1)
  if (nrow(row) == 0) {
    row <- tibble(term = "groupdiabetes_DN", estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_)
  }
  row
}

# Nest by protein and fit
results <- df_long %>%
  group_by(Protein) %>%
  group_nest(.key = "data") %>%
  mutate(tidy_row = map(data, fit_one)) %>%
  select(-data) %>%
  unnest(tidy_row) %>%
  transmute(
    Protein,
    estimate  = estimate,
    std_error = std.error,
    statistic = statistic,
    pvalue    = p.value
  ) %>%
  mutate(Padj = p.adjust(pvalue, method = "BH"))

# Group means for context
group_means <- df_long %>%
  group_by(Protein, group) %>%
  summarise(mean_NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_NPX, names_prefix = "mean_")

# Merge and order
final_results <- results %>%
  left_join(group_means, by = "Protein") %>%
  arrange(Padj)

# Save results
output_file <- file.path(out_dir, "tables", "protein_linear_regression_results.csv")
write.csv(final_results, output_file, row.names = FALSE)

message("Saved full regression results to: ", output_file)
print(head(final_results, 10))

# Count significant ones Padj<0.05
n_signif <- sum(final_results$Padj < 0.05, na.rm = TRUE)
cat("\nNumber of proteins with Padj < 0.05: ", n_signif
    , " out of ", nrow(final_results), "\n", sep = "")
cat("====================================================\n")



## -----------------------------------------------------------------------------
## NPX value range analysis
## -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

# Convert to numeric matrix (excluding eid)
olink_values <- olink_raw[, -1]
numeric_cols <- sapply(olink_values, is.numeric)
olink_values <- olink_values[, numeric_cols, drop = FALSE]

# Compute global range
olink_all <- unlist(olink_values)
olink_all <- olink_all[!is.na(olink_all)]
olink_min <- min(olink_all, na.rm = TRUE)
olink_max <- max(olink_all, na.rm = TRUE)
cat("OLINK NPX global range: min =", round(olink_min, 3), ", max =", round(olink_max, 3), "\n")

# Compute per-protein min and max
protein_ranges <- data.frame(
  Protein = colnames(olink_values),
  Min = apply(olink_values, 2, function(x) min(x, na.rm = TRUE)),
  Max = apply(olink_values, 2, function(x) max(x, na.rm = TRUE))
) %>%
  mutate(Range = Max - Min)

# Save summary table
range_file <- file.path(out_dir, "tables", "protein_NPX_ranges.csv")
write.csv(protein_ranges, range_file, row.names = FALSE)
message("Saved NPX range table to: ", range_file)

# Plot distributions of min and max
p_minmax <- protein_ranges %>%
  tidyr::pivot_longer(cols = c("Min", "Max"), names_to = "Type", values_to = "Value") %>%
  ggplot(aes(x = Value, fill = Type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("skyblue", "tomato")) +
  labs(
    title = "Distribution of NPX min and max values per protein",
    x = "NPX value",
    y = "Number of proteins",
    fill = ""
  ) +
  theme_bw(base_size = 13)

fig_file <- file.path(out_dir, "figs", "NPX_min_max_distribution.png")
ggsave(fig_file, p_minmax, width = 7, height = 5, dpi = 300)
message("Saved NPX min/max distribution plot to: ", fig_file)



# 
# ## -----------------------------------------------------------------------------
# ## Volcano plot: log2 fold change (DN vs nonDN) vs significance  [NA-filtered]
# ## -----------------------------------------------------------------------------
# library(ggplot2)
# 
# volcano_df <- final_results %>%
#   dplyr::filter(!is.na(Padj), !is.na(estimate)) %>%
#   dplyr::mutate(
#     Padj_plot     = pmax(Padj, 1e-300),               # avoid Inf in -log10
#     neg_log10Padj = -log10(Padj_plot),
#     sig           = ifelse(Padj < 0.05, "Significant", "Not significant")
#   )
# 
# p_volcano <- ggplot(volcano_df, aes(x = estimate, y = neg_log10Padj)) +
#   geom_point(aes(color = sig), alpha = 0.7, size = 2) +
#   scale_color_manual(values = c("gray60", "red")) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
#   geom_vline(xintercept = 0,            linetype = "solid",  color = "black") +
#   coord_cartesian(xlim = c(-2, 2)) +
#   labs(
#     title = "Volcano plot: diabetes_DN vs diabetes_nonDN",
#     x = "log2 fold change (DN / nonDN)",
#     y = "-log10(FDR-adjusted p-value)",
#     color = NULL
#   ) +
#   theme_bw(base_size = 13) +
#   theme(legend.position = "top")
# 
# fig_file <- file.path(out_dir, "figs", "volcano_diabetesDN_vs_nonDN.png")
# ggsave(fig_file, p_volcano, width = 7, height = 6, dpi = 300)
# message("Saved volcano plot to: ", fig_file)
# 
# 
# 
# ## -----------------------------------------------------------------------------
# ## Volcano plot with top 10 proteins labeled by absolute log2 fold-change [NA-filtered]
# ## -----------------------------------------------------------------------------
# library(ggplot2)
# library(ggrepel)
# 
# volcano_df <- final_results %>%
#   dplyr::filter(!is.na(Padj), !is.na(estimate)) %>%
#   dplyr::mutate(
#     Padj_plot     = pmax(Padj, 1e-300),
#     neg_log10Padj = -log10(Padj_plot),
#     sig           = ifelse(Padj < 0.05, "Significant", "Not significant")
#   )
# 
# top10_labels <- volcano_df %>%
#   dplyr::arrange(dplyr::desc(abs(estimate))) %>%
#   dplyr::slice_head(n = 10)
# 
# p_volcano_labeled <- ggplot(volcano_df, aes(x = estimate, y = neg_log10Padj)) +
#   geom_point(aes(color = sig), alpha = 0.7, size = 2) +
#   geom_text_repel(
#     data = top10_labels,
#     aes(label = Protein),
#     size = 3.5,
#     max.overlaps = 20,
#     box.padding = 0.4,
#     point.padding = 0.3,
#     segment.color = "gray50"
#   ) +
#   scale_color_manual(values = c("gray60", "red")) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
#   geom_vline(xintercept = 0,            linetype = "solid",  color = "black") +
#   coord_cartesian(xlim = c(-2, 2)) +
#   labs(
#     title = "Volcano plot: diabetes_DN vs diabetes_nonDN (top 10 |log2FC| labeled)",
#     x = "log2 fold change (DN / nonDN)",
#     y = "-log10(FDR-adjusted p-value)",
#     color = NULL
#   ) +
#   theme_bw(base_size = 13) +
#   theme(legend.position = "top")
# 
# fig_file_labeled <- file.path(out_dir, "figs", "volcano_diabetesDN_vs_nonDN_top10.png")
# ggsave(fig_file_labeled, p_volcano_labeled, width = 8, height = 6, dpi = 300)
# message("Saved volcano plot with top 10 labels to: ", fig_file_labeled)
# 
# 
# 
# ## -----------------------------------------------------------------------------
# ## Violin plots for top 10 proteins by |log2FC| (DN vs nonDN) in a 5x2 panel
# ## -----------------------------------------------------------------------------
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # 1) Decide the top 10 proteins by |estimate|
# top10_tbl <- final_results %>%
#   filter(!is.na(estimate), !is.na(Padj)) %>%
#   arrange(desc(abs(estimate))) %>%
#   slice_head(n = 10) %>%
#   select(Protein)
# 
# # 2) Build a long-form table of NPX values for those proteins
# #    If you already have df_long from earlier, this will be quick.
# if (!exists("df_long")) {
#   # Recreate df_long from combined; proteins start at column 28
#   meta_keep   <- c("record_id", "eid", "diabetes", "MNSI", "pn", "sex", "age_at_recruitment")
#   protein_cols <- colnames(combined)[28:ncol(combined)]
#   
#   analysis_df <- combined %>%
#     select(any_of(meta_keep), all_of(protein_cols)) %>%
#     filter(!is.na(diabetes), !is.na(MNSI), !is.na(pn), !is.na(sex), !is.na(age_at_recruitment)) %>%
#     mutate(group = case_when(
#       diabetes == 1 & pn == 1 ~ "diabetes_DN",
#       diabetes == 1 & pn == 0 ~ "diabetes_nonDN",
#       TRUE ~ NA_character_
#     )) %>%
#     filter(!is.na(group)) %>%
#     mutate(group = factor(group, levels = c("diabetes_nonDN", "diabetes_DN")))
#   
#   df_long <- analysis_df %>%
#     pivot_longer(cols = all_of(protein_cols), names_to = "Protein", values_to = "NPX")
# }
# 
# plot_df <- df_long %>%
#   semi_join(top10_tbl, by = "Protein") %>%
#   mutate(group = factor(group, levels = c("diabetes_nonDN", "diabetes_DN")))
# 
# # 3) Violin + boxplot + light jitter, faceted 5x2
# p_top10 <- ggplot(plot_df, aes(x = group, y = NPX)) +
#   geom_violin(aes(fill = group), trim = FALSE, scale = "width", alpha = 0.7) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
#   geom_jitter(width = 0.12, alpha = 0.15, size = 0.6) +
#   facet_wrap(~ Protein, ncol = 2, scales = "fixed") +
#   scale_fill_manual(values = c("diabetes_nonDN" = "gray70", "diabetes_DN" = "tomato")) +
#   labs(
#     title = "Top 10 proteins by |log2FC|: DN vs nonDN",
#     x = NULL,
#     y = "NPX (log2)"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     legend.position = "none",
#     strip.background = element_rect(fill = "white"),
#     axis.text.x = element_text(angle = 0, hjust = 0.5)
#   )
# 
# # 4) Save: 5 rows x 2 cols layout fits well at ~8 x 12 inches
# fig_file_violin <- file.path(out_dir, "figs", "top10_violin_DN_vs_nonDN_5x2.png")
# ggsave(fig_file_violin, p_top10, width = 8, height = 12, dpi = 300)
# message("Saved 5x2 violin panel to: ", fig_file_violin)
# 



# 
# ## =============================================================================
# ## Gene-set overlap and highlighting
# ## =============================================================================
# library(dplyr)
# library(readr)
# library(ggplot2)
# library(ggrepel)
# library(VennDetail)  # already loaded earlier via pacman, but load explicitly here
# 
# marker_dir <- file.path("rawdata", "markerSets")
# marker_files <- c(
#   "Stephanie-DRG.csv",
#   "Stephanie-Nerve.csv",
#   "Stephanie-SC.csv",
#   "Chen-Average-ORs.csv"
# )
# 
# ## -----------------------------------------------------------------------------
# ## Helpers
# ## -----------------------------------------------------------------------------
# read_genes_csv <- function(path, col = "Human") {
#   stopifnot(file.exists(path))
#   tb <- suppressMessages(read.csv(path, check.names = FALSE))
#   if (!col %in% names(tb)) stop("Column '", col, "' not found in ", path)
#   # Keep as character, trim spaces
#   genes <- unique(trimws(as.character(tb[[col]])))
#   genes[genes != "" & !is.na(genes)]
# }
# 
# # Build a clean volcano dataframe from final_results
# make_volcano_df <- function(final_results_df) {
#   final_results_df %>%
#     filter(!is.na(Padj), !is.na(estimate)) %>%
#     mutate(
#       Padj_plot     = pmax(Padj, 1e-300),
#       neg_log10Padj = -log10(Padj_plot),
#       sig           = ifelse(Padj < 0.05, "Significant", "Not significant")
#     )
# }
# 
# # Generic VennDetail plotting with graceful fallback
# plot_venn_upset <- function(named_lists, prefix, out_dir_figs) {
#   dir.create(out_dir_figs, showWarnings = FALSE, recursive = TRUE)
#   
#   # Build VennDetail object
#   ven <- VennDetail::venndetail(named_lists)
#   
#   # Traditional Venn
#   p_venn <- plot(ven, type = "venn")
#   venn_png <- file.path(out_dir_figs, paste0(prefix, "_venn.png"))
#   ggsave(venn_png, p_venn, width = 8, height = 6, dpi = 300)
#   message("Saved: ", venn_png)
#   
#   # UpSet
#   p_upset <- upset_plot(ven)
#   upset_png <- file.path(out_dir_figs, paste0(prefix, "_upset.png"))
#   ggsave(upset_png, p_upset, width = 8, height = 6, dpi = 300)
#   message("Saved: ", upset_png)
# }
# 
# # Volcano with a highlighted gene set
# plot_highlighted_volcano <- function(volcano_df, highlight_genes, title_tag, out_png, xlim = c(-2, 2)) {
#   highlight_genes <- unique(highlight_genes)
#   volc <- volcano_df %>%
#     mutate(
#       set_flag = ifelse(Protein %in% highlight_genes, "In set", "Other")
#     )
#   # Draw
#   p <- ggplot(volc, aes(x = estimate, y = neg_log10Padj)) +
#     geom_point(data = subset(volc, set_flag == "Other"),
#                color = "gray70", alpha = 0.6, size = 2) +
#     geom_point(data = subset(volc, set_flag == "In set"),
#                color = "dodgerblue3", alpha = 0.9, size = 2.2) +
#     geom_text_repel(
#       data = subset(volc, set_flag == "In set"),
#       aes(label = Protein),
#       size = 3.3,
#       max.overlaps = 50,
#       box.padding = 0.35,
#       point.padding = 0.25,
#       segment.color = "gray50"
#     ) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
#     geom_vline(xintercept = 0, linetype = "solid", color = "black") +
#     coord_cartesian(xlim = xlim) +
#     labs(
#       title = paste0("Volcano: DN vs nonDN  [", title_tag, " highlighted]"),
#       x = "log2 fold change (DN / nonDN)",
#       y = "-log10(FDR-adjusted p-value)"
#     ) +
#     theme_bw(base_size = 13) +
#     theme(legend.position = "none")
#   ggsave(out_png, p, width = 8, height = 6, dpi = 300)
#   message("Saved: ", out_png)
# }
# 
# ## -----------------------------------------------------------------------------
# ## 1) Load the four gene sets
# ## -----------------------------------------------------------------------------
# gene_sets <- list()
# for (f in marker_files) {
#   fpath <- file.path(marker_dir, f)
#   if (!file.exists(fpath)) stop("Gene set not found: ", fpath)
#   genes <- read_genes_csv(fpath, col = "Human")
#   key <- sub("\\.csv$", "", f)
#   gene_sets[[key]] <- genes
# }
# # Sanity
# # str(gene_sets)
# 
# ## -----------------------------------------------------------------------------
# ## 2) Overlap: significant list vs each set -> Venn and UpSet
# ## -----------------------------------------------------------------------------
# volcano_df <- make_volcano_df(final_results)
# sig_genes  <- volcano_df %>% filter(Padj < 0.05) %>% pull(Protein) %>% unique()
# 
# # Three Stephanie sets
# steph_sets <- list(
#   Stephanie_DRG   = gene_sets[["Stephanie-DRG"]],
#   Stephanie_Nerve = gene_sets[["Stephanie-Nerve"]],
#   Stephanie_SC    = gene_sets[["Stephanie-SC"]]
# )
# 
# # IEU-OpenGWAS set from Chen
# ieu_openGWAS <- gene_sets[["Chen-Average-ORs"]]
# 
# # Venn + UpSet for significant vs Stephanie sets (combine as one list)
# named_lists_steph <- c(list(Significant = sig_genes), steph_sets)
# plot_venn_upset(named_lists_steph, prefix = "sig_vs_Stephanie", out_dir_figs = file.path(out_dir, "figs"))
# 
# # Venn + UpSet for significant vs IEU-OpenGWAS
# named_lists_ieu <- list(Significant = sig_genes, IEU_OpenGWAS = ieu_openGWAS)
# plot_venn_upset(named_lists_ieu, prefix = "sig_vs_IEU-OpenGWAS", out_dir_figs = file.path(out_dir, "figs"))
# 
# ## -----------------------------------------------------------------------------
# ## 3) Senescence set = union of the three Stephanie files
# ##    Highlight and label all senescence proteins that exist in results
# ## -----------------------------------------------------------------------------
# senescence_genes <- unique(unlist(steph_sets))
# senescence_in_results <- intersect(senescence_genes, volcano_df$Protein)
# 
# png_senescence <- file.path(out_dir, "figs", "volcano_highlight_senescence.png")
# plot_highlighted_volcano(volcano_df, senescence_in_results, "Senescence (Stephanie union)", png_senescence)
# 
# ## -----------------------------------------------------------------------------
# ## 4) IEU-OpenGWAS highlighting as a separate volcano
# ## -----------------------------------------------------------------------------
# ieu_in_results <- intersect(ieu_openGWAS, volcano_df$Protein)
# png_ieu <- file.path(out_dir, "figs", "volcano_highlight_IEU-OpenGWAS.png")
# plot_highlighted_volcano(volcano_df, ieu_in_results, "IEU-OpenGWAS", png_ieu)
# 
# ## -----------------------------------------------------------------------------
# ## 5) Save overlap tables for record  
# ## -----------------------------------------------------------------------------
# overlap_dir <- file.path(out_dir, "tables")
# dir.create(overlap_dir, showWarnings = FALSE, recursive = TRUE)
# 
# write_overlap <- function(gene_vec, sig_vec, out_path) {
#   df <- data.frame(Protein = sort(unique(gene_vec)), stringsAsFactors = FALSE)
#   df$In_Significant <- df$Protein %in% sig_vec
#   write.csv(df, out_path, row.names = FALSE)
#   message("Saved table: ", out_path)
# }
# 
# write_overlap(senescence_in_results, sig_genes,
#               file.path(overlap_dir, "senescence_overlap_in_results.csv"))
# 
# write_overlap(ieu_in_results, sig_genes,
#               file.path(overlap_dir, "IEU-OpenGWAS_overlap_in_results.csv"))
# 
# 
# 



## =============================================================================
## Gene-set overlap and highlighting - regular Venn and UpSetR
## =============================================================================
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(UpSetR)
library(tidyr)

marker_dir <- file.path("rawdata", "markerSets")
marker_files <- c(
  "Stephanie-DRG.csv",
  "Stephanie-Nerve.csv",
  "Stephanie-SC.csv",
  "Chen-Average-ORs.csv",
  "HumanGeneSymbol-MouseGeneSymbol.csv"  # SenMayo panel
)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
read_genes_csv <- function(path, col = "Human") {
  stopifnot(file.exists(path))
  tb <- suppressMessages(read.csv(path, check.names = FALSE))
  if (!col %in% names(tb)) stop("Column '", col, "' not found in ", path)
  genes <- unique(trimws(as.character(tb[[col]])))
  genes[genes != "" & !is.na(genes)]
}

make_volcano_df <- function(final_results_df) {
  final_results_df %>%
    filter(!is.na(Padj), !is.na(estimate)) %>%
    mutate(
      Padj_plot     = pmax(Padj, 1e-300),
      neg_log10Padj = -log10(Padj_plot),
      sig           = ifelse(Padj < 0.05, "Significant", "Not significant")
    )
}

# Simple ggVennDiagram with custom colors
plot_venn_regular <- function(named_lists, prefix, out_dir_figs, set_colors = NULL) {
  dir.create(out_dir_figs, showWarnings = FALSE, recursive = TRUE)
  # default colors if not provided
  if (is.null(set_colors)) {
    set_colors <- setNames(
      c("#1f77b4","#ff7f0e","#2ca02c","#9467bd","#8c564b")[seq_along(named_lists)],
      names(named_lists)
    )
  }
  p <- ggVennDiagram(
    named_lists,
    label = "count",
    label_alpha = 0,
    set_color = set_colors
  ) + scale_fill_gradient(low = "white", high = "red")
  out <- file.path(out_dir_figs, paste0(prefix, "_ggVennDiagram.png"))
  ggsave(out, p, width = 8, height = 6, dpi = 300)
  message("Saved Venn: ", out)
}

# UpSetR plot from a named list
plot_upset_regular <- function(named_lists, prefix, out_dir_figs, set_colors = NULL, nintersects = 40) {
  dir.create(out_dir_figs, showWarnings = FALSE, recursive = TRUE)
  if (is.null(set_colors)) {
    set_colors <- setNames(
      c("#1f77b4","#ff7f0e","#2ca02c","#9467bd","#8c564b")[seq_along(named_lists)],
      names(named_lists)
    )
  }
  df_up <- UpSetR::fromList(named_lists)
  # ensure the sets bar colors align with column order
  ordered_colors <- set_colors[colnames(df_up)]
  # UpSetR uses base graphics; save via png device
  png(filename = file.path(out_dir_figs, paste0(prefix, "_UpSetR.png")),
      width = 1800, height = 800, res = 150)
  print(
    UpSetR::upset(
      df_up,
      nsets = ncol(df_up),
      nintersects = nintersects,
      sets = colnames(df_up),
      sets.x.label = "Set size",
      mainbar.y.label = "Intersection size",
      point.size = 3,
      line.size = 0.8,
      sets.bar.color = unname(ordered_colors),
      order.by = "freq"
    )
  )
  dev.off()
  message("Saved UpSet: ", file.path(out_dir_figs, paste0(prefix, "_UpSetR.png")))
}

# # Volcano highlighting helper
# plot_highlighted_volcano <- function(volcano_df, highlight_genes, title_tag, out_png, xlim = c(-2, 2)) {
#   highlight_genes <- unique(highlight_genes)
#   volc <- volcano_df %>%
#     mutate(set_flag = ifelse(Protein %in% highlight_genes, "In set", "Other"))
#   p <- ggplot(volc, aes(x = estimate, y = neg_log10Padj)) +
#     geom_point(data = subset(volc, set_flag == "Other"),
#                color = "gray70", alpha = 0.6, size = 2) +
#     geom_point(data = subset(volc, set_flag == "In set"),
#                color = "dodgerblue3", alpha = 0.9, size = 2.2) +
#     geom_text_repel(
#       data = subset(volc, set_flag == "In set"),
#       aes(label = Protein),
#       size = 3.3,
#       max.overlaps = 50,
#       box.padding = 0.35,
#       point.padding = 0.25,
#       segment.color = "gray50"
#     ) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
#     geom_vline(xintercept = 0, linetype = "solid", color = "black") +
#     coord_cartesian(xlim = xlim) +
#     labs(
#       title = paste0("Volcano: DN vs nonDN [", title_tag, " highlighted]"),
#       x = "log2 fold change (DN / nonDN)",
#       y = "-log10(FDR-adjusted p-value)"
#     ) +
#     theme_bw(base_size = 13) +
#     theme(legend.position = "none")
#   ggsave(out_png, p, width = 8, height = 6, dpi = 300)
#   message("Saved: ", out_png)
# }

# Volcano highlighting helper (only highlight significant ones)
plot_highlighted_volcano <- function(volcano_df, highlight_genes, title_tag, out_png, xlim = c(-2, 2)) {
  highlight_genes <- unique(highlight_genes)
  
  volc <- volcano_df %>%
    mutate(
      sig_flag = Padj < 0.05,
      set_flag = case_when(
        Protein %in% highlight_genes & sig_flag ~ "In set (sig)",
        TRUE ~ "Other"
      )
    )
  
  p <- ggplot(volc, aes(x = estimate, y = neg_log10Padj)) +
    geom_point(data = subset(volc, set_flag == "Other"),
               color = "gray70", alpha = 0.6, size = 2) +
    geom_point(data = subset(volc, set_flag == "In set (sig)"),
               color = "dodgerblue3", alpha = 0.9, size = 2.2) +
    geom_text_repel(
      data = subset(volc, set_flag == "In set (sig)"),
      aes(label = Protein),
      size = 3.3,
      max.overlaps = 50,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "gray50"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    coord_cartesian(xlim = xlim) +
    labs(
      title = paste0("Volcano: DN vs nonDN [", title_tag, " highlighted]"),
      x = "log2 fold change (DN / nonDN)",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none")
  
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  message("Saved: ", out_png)
}


# Simple writer for overlap tables
write_overlap <- function(gene_vec, sig_vec, out_path) {
  df <- data.frame(Protein = sort(unique(gene_vec)), stringsAsFactors = FALSE)
  df$In_Significant <- df$Protein %in% sig_vec
  write.csv(df, out_path, row.names = FALSE)
  message("Saved table: ", out_path)
}

## -----------------------------------------------------------------------------
## 1) Load the four gene sets
## -----------------------------------------------------------------------------
gene_sets <- list()
for (f in marker_files) {
  fpath <- file.path(marker_dir, f)
  if (!file.exists(fpath)) stop("Gene set not found: ", fpath)
  genes <- read_genes_csv(fpath, col = "Human")
  key <- sub("\\.csv$", "", f)
  gene_sets[[key]] <- genes
}

## Change the name for HumanGeneSymbol-MouseGeneSymbol" -> "SenMayo"
names(gene_sets)[names(gene_sets) == "HumanGeneSymbol-MouseGeneSymbol"] <- "SenMayo"


## -----------------------------------------------------------------------------
## 2) Build significant list from results
## -----------------------------------------------------------------------------
volcano_df <- make_volcano_df(final_results)
sig_genes  <- volcano_df %>% filter(Padj < 0.05) %>% pull(Protein) %>% unique()

## -----------------------------------------------------------------------------
## 3) Senescence set = union of Stephanie files; IEU set from Chen
## -----------------------------------------------------------------------------
steph_sets <- list(
  Stephanie_DRG   = gene_sets[["Stephanie-DRG"]],
  Stephanie_Nerve = gene_sets[["Stephanie-Nerve"]],
  Stephanie_SC    = gene_sets[["Stephanie-SC"]]
)
senescence_genes   <- unique(unlist(steph_sets))
ieu_openGWAS_genes <- gene_sets[["Chen-Average-ORs"]]
senmayo_genes   <- gene_sets[["SenMayo"]]


senescence_in_results <- intersect(senescence_genes, volcano_df$Protein)
ieu_in_results        <- intersect(ieu_openGWAS_genes, volcano_df$Protein)
senmayo_in_results        <- intersect(senmayo_genes, volcano_df$Protein)

## -----------------------------------------------------------------------------
## 4) Venn and UpSet for overlaps
## -----------------------------------------------------------------------------
venn_colors <- c(
  Significant = "#333333",
  Stephanie_DRG = "#1f77b4",
  Stephanie_Nerve = "#2ca02c",
  Stephanie_SC = "#9467bd",
  IEU_OpenGWAS = "#ff7f0e"
)

# Significant vs Stephanie sets
named_lists_steph <- c(list(Significant = sig_genes), steph_sets)
plot_venn_regular(named_lists_steph, prefix = "sig_vs_Stephanie", out_dir_figs = file.path(out_dir, "figs"),
                  set_colors = venn_colors[names(named_lists_steph)])
plot_upset_regular(named_lists_steph, prefix = "sig_vs_Stephanie", out_dir_figs = file.path(out_dir, "figs"),
                   set_colors = venn_colors[names(named_lists_steph)], nintersects = 40)

# Significant vs IEU
named_lists_ieu <- list(Significant = sig_genes, IEU_OpenGWAS = ieu_openGWAS_genes)
plot_venn_regular(named_lists_ieu, prefix = "sig_vs_IEU-OpenGWAS", out_dir_figs = file.path(out_dir, "figs"),
                  set_colors = venn_colors[names(named_lists_ieu)])
plot_upset_regular(named_lists_ieu, prefix = "sig_vs_IEU-OpenGWAS", out_dir_figs = file.path(out_dir, "figs"),
                   set_colors = venn_colors[names(named_lists_ieu)], nintersects = 40)

## -----------------------------------------------------------------------------
## 5) Highlighted volcano plots for senescence and IEU sets
## -----------------------------------------------------------------------------
png_senescence <- file.path(out_dir, "figs", "volcano_highlight_senescence.png")
plot_highlighted_volcano(volcano_df, senescence_in_results, "Senescence (Stephanie union)", png_senescence)

png_ieu <- file.path(out_dir, "figs", "volcano_highlight_IEU-OpenGWAS.png")
plot_highlighted_volcano(volcano_df, ieu_in_results, "IEU-OpenGWAS", png_ieu)

png_senmayo <- file.path(out_dir, "figs", "volcano_highlight_SenMayo.png")
plot_highlighted_volcano(volcano_df, senmayo_in_results, "SenMayo", png_senmayo)

## -----------------------------------------------------------------------------
## 6) Save overlap tables
## -----------------------------------------------------------------------------
overlap_dir <- file.path(out_dir, "tables")
dir.create(overlap_dir, showWarnings = FALSE, recursive = TRUE)

write_overlap(senescence_in_results, sig_genes,
              file.path(overlap_dir, "senescence_overlap_in_results.csv"))

write_overlap(ieu_in_results, sig_genes,
              file.path(overlap_dir, "IEU-OpenGWAS_overlap_in_results.csv"))

write_overlap(senmayo_in_results, sig_genes,
              file.path(overlap_dir, "SenMayo_overlap_in_results.csv"))





## =============================================================================
## Violin panels for overlapping genes: Senescence and IEU-OpenGWAS
## =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Helper to ensure df_long exists (proteins start at column 28 in `combined`)
ensure_df_long <- function() {
  if (exists("df_long")) return(invisible(NULL))
  meta_keep   <- c("record_id", "eid", "diabetes", "MNSI", "pn", "sex", "age_at_recruitment")
  protein_cols <- colnames(combined)[28:ncol(combined)]
  
  analysis_df <<- combined %>%
    select(any_of(meta_keep), all_of(protein_cols)) %>%
    filter(!is.na(diabetes), !is.na(MNSI), !is.na(pn), !is.na(sex), !is.na(age_at_recruitment)) %>%
    mutate(group = case_when(
      diabetes == 1 & pn == 1 ~ "diabetes_DN",
      diabetes == 1 & pn == 0 ~ "diabetes_nonDN",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group)) %>%
    mutate(group = factor(group, levels = c("diabetes_nonDN", "diabetes_DN")))
  
  df_long <<- analysis_df %>%
    pivot_longer(cols = all_of(protein_cols), names_to = "Protein", values_to = "NPX")
}

# Helper to order facets by significance then effect size
facet_order_from_results <- function(genes_vec) {
  volcano_df %>%
    filter(Protein %in% genes_vec) %>%
    arrange(Padj, desc(abs(estimate))) %>%
    pull(Protein) %>%
    unique()
}

# Generic violin panel maker
plot_violin_for_genes <- function(genes_vec, tag, file_out, ncol = 3) {
  genes_vec <- unique(genes_vec)
  if (length(genes_vec) == 0) {
    message("No overlapping genes for ", tag, ". Skipping violin plot.")
    return(invisible(NULL))
  }
  ensure_df_long()
  ord <- facet_order_from_results(genes_vec)
  plot_df <- df_long %>%
    filter(Protein %in% genes_vec) %>%
    mutate(
      group = factor(group, levels = c("diabetes_nonDN", "diabetes_DN")),
      Protein = factor(Protein, levels = ord)
    )
  
  # --- Plot ---
  p <- ggplot(plot_df, aes(x = group, y = NPX)) +
    geom_violin(aes(fill = group), trim = FALSE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.12, alpha = 0.15, size = 0.6) +
    facet_wrap(~ Protein, ncol = ncol, scales = "fixed") +
    scale_fill_manual(values = c("diabetes_nonDN" = "gray70", "diabetes_DN" = "tomato")) +
    labs(
      title = paste0(tag, " overlap: DN vs nonDN"),
      x = NULL,
      y = "NPX (log2)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  # --- Size control ---
  height_in <- min(max(6, ceiling(length(genes_vec) / ncol) * 2.5), 50)  # cap at 50 inches
  width_in  <- min(10, 50)
  
  ggsave(
    file_out,
    p,
    width = width_in,
    height = height_in,
    dpi = 300,
    limitsize = FALSE  # allow large plots safely
  )
  message("Saved violin panel: ", file_out, " (", length(genes_vec), " genes)")
}


# Build overlap vectors from earlier steps
# senescence_in_results and ieu_in_results were defined in the previous block.
# If not, reconstruct them quickly:
if (!exists("senescence_in_results") || !exists("ieu_in_results")) {
  # Requires steph_sets, gene_sets, and volcano_df from the previous block
  if (!exists("steph_sets")) {
    steph_sets <- list(
      Stephanie_DRG   = gene_sets[["Stephanie-DRG"]],
      Stephanie_Nerve = gene_sets[["Stephanie-Nerve"]],
      Stephanie_SC    = gene_sets[["Stephanie-SC"]]
    )
  }
  if (!exists("volcano_df")) {
    volcano_df <- final_results %>%
      filter(!is.na(Padj), !is.na(estimate)) %>%
      mutate(Padj_plot = pmax(Padj, 1e-300), neg_log10Padj = -log10(Padj_plot))
  }
  senescence_genes <- unique(unlist(steph_sets))
  senescence_in_results <- intersect(senescence_genes, volcano_df$Protein)
  ieu_openGWAS_genes <- gene_sets[["Chen-Average-ORs"]]
  ieu_in_results <- intersect(ieu_openGWAS_genes, volcano_df$Protein)
}

# Output files
fig_senescence <- file.path(out_dir, "figs", "violin_senescence_overlap_DN_vs_nonDN.png")
fig_ieu        <- file.path(out_dir, "figs", "violin_IEU-OpenGWAS_overlap_DN_vs_nonDN.png")

# Create the two panels
plot_violin_for_genes(senescence_in_results, "Senescence (Stephanie union)", fig_senescence, ncol = 2)
plot_violin_for_genes(ieu_in_results, "IEU-OpenGWAS", fig_ieu, ncol = 2)


# If not, also reconstruct SenMayo quickly:
if (!exists("senmayo_in_results")) {
  if (!exists("gene_sets")) stop("gene_sets object not found for SenMayo reconstruction")
  if (!exists("volcano_df")) {
    volcano_df <- final_results %>%
      filter(!is.na(Padj), !is.na(estimate)) %>%
      mutate(Padj_plot = pmax(Padj, 1e-300), neg_log10Padj = -log10(Padj_plot))
  }
  senmayo_genes <- gene_sets[["SenMayo"]]
  senmayo_in_results <- intersect(senmayo_genes, volcano_df$Protein)
}

# Output file for SenMayo
fig_senmayo <- file.path(out_dir, "figs", "violin_SenMayo_overlap_DN_vs_nonDN.png")

# Create the SenMayo panel
plot_violin_for_genes(senmayo_in_results, "SenMayo", fig_senmayo, ncol = 2)










## Save the current session data with compression
save.image(file = file.path(out_dir, "plasma_proteomics_analysis_workspace.RData"), compress = "gzip")



## -----------------------------------------------------------------------------
## End of script

