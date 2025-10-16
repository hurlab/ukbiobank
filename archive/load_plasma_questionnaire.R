#!/usr/bin/env Rscript

# Combine OLINK protein measurements with MNSI questionnaire responses.
# MNSI scoring is performed inline. If the raw MNSI is NA due to missing items,
# we keep the partial sum when the sum of available (non-NA) scored items is ≥ 3.

options(stringsAsFactors = FALSE)

olink_path    <- file.path("rawdata", "plasma", "OLINK protein with name_instance0.csv")
question_path <- file.path("rawdata", "plasma", "Questionnaire 120071 to 120085.csv")
output_path   <- file.path("rawdata", "plasma", "olink_questionnaire_matched.csv")

# ---------- I/O checks ----------
if (!file.exists(olink_path)) {
  stop("OLINK protein file not found: ", olink_path)
}
if (!file.exists(question_path)) {
  stop("Questionnaire file not found: ", question_path)
}

# ---------- Load OLINK ----------
olink_raw <- read.csv(olink_path, check.names = FALSE)
colnames(olink_raw)[1] <- "eid"
protein_labels  <- colnames(olink_raw)
protein_symbols <- protein_labels
protein_symbols[-1] <- trimws(sub(";.*$", "", protein_labels[-1]))
colnames(olink_raw) <- make.unique(protein_symbols)
olink_raw$eid <- as.character(olink_raw$eid)

# ---------- Load questionnaire ----------
question_raw <- read.csv(question_path, check.names = FALSE)
colnames(question_raw)[1] <- "record_id"
colnames(question_raw)[ncol(question_raw)] <- "eid"
question_codes <- 120071:120085
question_cols  <- sprintf("f.%06d.0.0", question_codes)
colnames(question_raw)[2:(length(question_cols) + 1)] <- question_cols
question_raw$eid <- as.character(question_raw$eid)

# Map responses to 0/1. Accept "YES"/"NO" and 0/1 already present.
map_response <- function(x) {
  # If numeric with 0/1, return as is
  if (is.numeric(x)) {
    return(ifelse(x %in% c(0, 1), as.numeric(x), NA_real_))
  }
  # Otherwise parse text
  nx <- trimws(toupper(as.character(x)))
  out <- rep(NA_real_, length(nx))
  out[nx == "YES"] <- 1
  out[nx == "NO"]  <- 0
  as.numeric(out)
}
question_raw[question_cols] <- lapply(question_raw[question_cols], map_response)

# Diabetes flag (UKB 120007). If missing, default to 1 as in your original script.
if (!"f.120007.0.0" %in% names(question_raw)) {
  question_raw$f.120007.0.0 <- 1L
}

# Drop duplicate questionnaire rows by eid, keep first
dup_idx <- duplicated(question_raw$eid)
if (any(dup_idx)) {
  warning("Duplicate questionnaire IDs detected; keeping first occurrence.")
  question_raw <- question_raw[!dup_idx, ]
}

# ---------- Match participants present in both datasets ----------
common_ids <- sort(intersect(question_raw$eid, olink_raw$eid))
if (!length(common_ids)) {
  stop("No participants found in both datasets.")
}
question_matched <- question_raw[match(common_ids, question_raw$eid), ]
olink_matched    <- olink_raw[match(common_ids, olink_raw$eid), ]

missing_rows <- is.na(question_matched$eid) | is.na(olink_matched$eid)
if (any(missing_rows)) {
  question_matched <- question_matched[!missing_rows, , drop = FALSE]
  olink_matched    <- olink_matched[!missing_rows, , drop = FALSE]
}

# ---------- MNSI scoring inside this script ----------
# Items 120071 to 120085 are the questionnaire. Two are reversed (120077, 120083).
# Two items are not scored in the standard self-admin MNSI: Q4 (120074) and Q10 (120080).
# We compute a raw sum across scored items. If that raw sum is NA only due to missingness,
# we keep the partial sum when the partial sum is ≥ 3.

# Ensure negatives are treated as missing
clean_negatives <- function(v) { v[v < 0] <- NA; v }

mnsi_all_cols <- sprintf("f.%06d.0.0", 120071:120085)
present_cols  <- intersect(mnsi_all_cols, names(question_matched))
question_matched[present_cols] <- lapply(question_matched[present_cols], clean_negatives)

# Build reversed items
if ("f.120083.0.0" %in% names(question_matched)) {
  question_matched$f.120083.0.0_rev <- ifelse(question_matched$f.120083.0.0 == 1, 0,
                                              ifelse(question_matched$f.120083.0.0 == 0, 1, NA))
} else {
  question_matched$f.120083.0.0_rev <- NA_real_
}
if ("f.120077.0.0" %in% names(question_matched)) {
  question_matched$f.120077.0.0_rev <- ifelse(question_matched$f.120077.0.0 == 1, 0,
                                              ifelse(question_matched$f.120077.0.0 == 0, 1, NA))
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

# Initialize outputs
question_matched$diabetes <- ifelse(!is.na(question_matched$f.120007.0.0) &
                                      question_matched$f.120007.0.0 == 1, 1L, 0L)
question_matched$MNSI <- NA_real_
question_matched$pn   <- NA_integer_

# Compute MNSI only for diabetes == 1, leave NA otherwise
idx <- which(question_matched$diabetes == 1)

if (length(idx)) {
  # Compact view of scored items
  M <- question_matched[idx, mnsi_items, drop = FALSE]
  
  # Sums and missingness
  raw_sum           <- rowSums(M, na.rm = FALSE)   # NA if any NA present
  partial_sum       <- rowSums(M, na.rm = TRUE)    # ignores NAs
  non_missing_count <- rowSums(!is.na(M))
  all_missing       <- non_missing_count == 0
  has_missing       <- non_missing_count < ncol(M)
  
  # Start with raw sums for complete cases
  final_sum <- raw_sum
  
  # For rows with any missing items:
  # - if partial_sum > 3, keep partial_sum
  # - else set to NA (so these rows are excluded downstream)
  keep_partial <- has_missing & (partial_sum > 3)
  drop_partial <- has_missing & !keep_partial
  
  final_sum[keep_partial] <- partial_sum[keep_partial]
  final_sum[drop_partial] <- NA_real_
  final_sum[all_missing]  <- NA_real_
  
  # Write back
  question_matched$MNSI[idx] <- final_sum
  
  # PN outcome: > 3. Rows with MNSI = NA remain NA here, so they do not fall into the ≤ 3 group.
  question_matched$pn[idx] <- ifelse(question_matched$MNSI[idx] > 3, 1L, 0L)
}


# ---------- Combine and write ----------
combined <- cbind(question_matched, olink_matched[, -1, drop = FALSE])
write.csv(combined, output_path, row.names = FALSE)
message("Wrote matched protein and questionnaire data to ", output_path)
