library(data.table)
library(dplyr)

# This script prepares a single symmetric geno_corr_df RData object from
# two inputs placed under ../inputs/:
#  - finngen_R13_*heritability.tsv (columns: PHENO, H2, SE, INT, INT_SE, RATIO, RATIO_SE)
#  - finngen_R13*.summary.tsv (ldsc rg summary with columns including p1, p2, rg, se, z, p,
#    h2_obs, h2_obs_se, h2_int, h2_int_se, gcov_int, gcov_int_se, CONVERGED)
#
# Output:
#  - ../Rdata_outputs/geno_correlation_sig.Rdata containing geno_corr_df keyed by (p1, p2)

inputs_dir <- "../inputs"

# Discover files
herit_file <- list.files(inputs_dir, pattern = "finngen_R13.*heritability\\.tsv$", full.names = TRUE)
summary_file <- list.files(inputs_dir, pattern = "finngen_R13.*summary\\.tsv$", full.names = TRUE)

if (length(herit_file) == 0) {
  stop("Could not find finngen_R13_*heritability.tsv under ../inputs")
}
if (length(summary_file) == 0) {
  stop("Could not find finngen_R13*.summary.tsv under ../inputs")
}

herit_file <- herit_file[1]
summary_file <- summary_file[1]

message("Using heritability file: ", herit_file)
message("Using rg summary file: ", summary_file)

# Load heritability table and normalize column names to expected h2_* fields for p1
h2_dt <- fread(herit_file)
setnames(h2_dt, old = c("PHENO", "H2", "SE", "INT", "INT_SE"),
               new = c("pheno", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se"),
               skip_absent = TRUE)

required_h2_cols <- c("pheno", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se")
missing_h2 <- setdiff(required_h2_cols, names(h2_dt))
if (length(missing_h2) > 0) {
  stop(paste0("Missing required columns in heritability tsv: ", paste(missing_h2, collapse = ", ")))
}

h2_map <- h2_dt[, .(pheno = as.character(pheno),
                    h2_obs = as.numeric(h2_obs),
                    h2_obs_se = as.numeric(h2_obs_se),
                    h2_int = as.numeric(h2_int),
                    h2_int_se = as.numeric(h2_int_se))]
setkey(h2_map, pheno)

# Load pairwise genetic correlation summary
rg_dt <- fread(summary_file)

# Ensure required columns exist
required_rg_cols <- c("p1", "p2", "rg", "se", "z", "p",
                      "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                      "gcov_int", "gcov_int_se")
missing_rg <- setdiff(required_rg_cols, names(rg_dt))
if (length(missing_rg) > 0) {
  stop(paste0("Missing required columns in summary tsv: ", paste(missing_rg, collapse = ", ")))
}

# Normalize types
rg_dt[, `:=`(
  p1 = as.character(p1),
  p2 = as.character(p2),
  rg = as.numeric(rg),
  se = as.numeric(se),
  z = as.numeric(z),
  p = as.numeric(p),
  h2_obs = as.numeric(h2_obs),
  h2_obs_se = as.numeric(h2_obs_se),
  h2_int = as.numeric(h2_int),
  h2_int_se = as.numeric(h2_int_se),
  gcov_int = as.numeric(gcov_int),
  gcov_int_se = as.numeric(gcov_int_se)
)]

# Add description columns if absent; default to IDs
if (!("description_p1" %in% names(rg_dt))) rg_dt[, description_p1 := p1]
if (!("description_p2" %in% names(rg_dt))) rg_dt[, description_p2 := p2]

## Build symmetric table for the Shiny browser (geno_corr_df)
rg_swapped <- copy(rg_dt)
rg_swapped[, `:=`(
  p1 = rg_dt$p2,
  p2 = rg_dt$p1,
  description_p1 = rg_dt$description_p2,
  description_p2 = rg_dt$description_p1
)]

# Join heritability values for p1 for both original and swapped tables to ensure correctness
rg_dt_joined <- merge(rg_dt[, -c("h2_obs", "h2_obs_se", "h2_int", "h2_int_se")],
                      h2_map, by.x = "p1", by.y = "pheno", all.x = TRUE)
rg_swapped_joined <- merge(rg_swapped[, -c("h2_obs", "h2_obs_se", "h2_int", "h2_int_se")],
                           h2_map, by.x = "p1", by.y = "pheno", all.x = TRUE)

geno_corr_df <- rbind(rg_dt_joined, rg_swapped_joined, fill = TRUE)
setkeyv(geno_corr_df, c("p1", "p2"))

# Save RData for Shiny (both in root and under r2_browser for app convenience)
dir.create("../Rdata_outputs", showWarnings = FALSE)
dir.create("../r2_browser/Rdata_outputs", showWarnings = FALSE)
save(geno_corr_df, file = "../Rdata_outputs/geno_correlation_sig.Rdata")
save(geno_corr_df, file = "../r2_browser/Rdata_outputs/geno_correlation_sig.Rdata")
message("Saved RData to ../Rdata_outputs and ../r2_browser/Rdata_outputs")

## Build lower-triangular .r2 file for static heatmaps
phens <- sort(unique(c(rg_dt$p1, rg_dt$p2)))

# Build description map per phenotype (prefer description_p1 if present)
desc_map <- unique(rbind(
  rg_dt[, .(pheno = p1, description = description_p1)],
  rg_dt[, .(pheno = p2, description = description_p2)]
))
setkey(desc_map, pheno)

rows <- vector("list", length = length(phens) * (length(phens) - 1L) / 2L)
idx <- 1L
for (i in seq_len(length(phens) - 1L)) {
  for (j in (i + 1L):length(phens)) {
    a <- phens[i]
    b <- phens[j]
    hit <- rg_dt[(p1 == a & p2 == b) | (p1 == b & p2 == a)]
    if (nrow(hit) < 1L) next
    hit <- hit[1]
    rows[[idx]] <- data.table(
      p1 = a,
      p2 = b,
      description_p1 = if (!is.na(desc_map[a, description])) desc_map[a, description] else a,
      description_p2 = if (!is.na(desc_map[b, description])) desc_map[b, description] else b,
      rg = as.numeric(hit$rg),
      se = as.numeric(hit$se),
      z = as.numeric(hit$z),
      p = as.numeric(hit$p),
      h2_obs = as.numeric(h2_map[b, h2_obs])
    )
    idx <- idx + 1L
  }
}

tri_dt <- rbindlist(rows)
dir.create("../r2_results", showWarnings = FALSE)
fwrite(tri_dt, file = "../r2_results/geno_correlation_sig.r2", sep = "\t", quote = FALSE, na = "NA")
message("Saved ../r2_results/geno_correlation_sig.r2 with ", nrow(tri_dt), " rows")

