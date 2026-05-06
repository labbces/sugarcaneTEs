# -------------------------------------------------------------
#                             Description
# R script to identify the 4 most abundant TEs by exon type
# and calculate the total length of the associated regions (pbCAZymes)
# Outuput file: "summary_Top4_TEs_with_size_per_species.tsv"
# -------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)

# -----------------------------
# SETTINGS
# -----------------------------

species <- c("c19", "c57", "cK3", "c85", "cNp", "c87", "cPur", "c80")
exon_types <- c("monoEXON", "firstEXON", "internalEXON", "finalEXON")

base_path <- "/home/R/CAZymeTEs_exonsType"
bad_classes <- c("Unknown", "Simple_repeat", "Low_complexity", "Satellite", ".")

# -----------------------------
# FUNCTION: top 4 TEs by count and total size
# -----------------------------

get_topTEs_size <- function(file_path) {
  if (!file.exists(file_path))
    return(data.frame(TE_type = rep(NA, 4), Count = rep(NA, 4), pbCAZymes = rep(NA, 4)))
  
  bed <- fread(file_path, sep = "\t", header = FALSE, fill = TRUE, quote = "")
  if (ncol(bed) < 16)
    return(data.frame(TE_type = rep(NA, 4), Count = rep(NA, 4), pbCAZymes = rep(NA, 4)))
  
  TE_type <- bed[[16]]
  start <- suppressWarnings(as.numeric(bed[[2]]))
  end <- suppressWarnings(as.numeric(bed[[3]]))
  
  # Sequence length (end - start + 1)
  len <- end - start + 1
  
  valid <- !(TE_type %in% bad_classes) & !is.na(TE_type) & !is.na(len) & len > 0
  bed_valid <- data.frame(TE_type = TE_type[valid], len = len[valid])
  
  if (nrow(bed_valid) == 0)
    return(data.frame(TE_type = rep(NA, 4), Count = rep(NA, 4), pbCAZymes = rep(NA, 4)))
  
  # Count and total length by TE type
  summary_tbl <- bed_valid %>%
    group_by(TE_type) %>%
    summarise(Count = n(), pbCAZymes = sum(len, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Count)) %>%
    head(4)
  
  # Ensure 4 rows even when fewer than 4 valid TE types are present
  if (nrow(summary_tbl) < 4) {
    summary_tbl <- rbind(summary_tbl,
                         data.frame(
                           TE_type = rep(NA, 4 - nrow(summary_tbl)),
                           Count = rep(NA, 4 - nrow(summary_tbl)),
                           pbCAZymes = rep(NA, 4 - nrow(summary_tbl))
                         ))
  }
  
  return(summary_tbl)
}

# -----------------------------
# MAIN LOOP
# -----------------------------

results_list <- list()

for (etype in exon_types) {
  message("Processing exon type: ", etype)
  
  # 4 rows per exon type (for the top 4 TEs)
  exon_result <- data.frame(ExonType = rep(etype, 4))
  
  for (sp in species) {
    fpath <- file.path(base_path, paste0(sp, "_", etype, ".closest_TEs-CAZymes.gtf"))
    topTEs <- get_topTEs_size(fpath)
    
    # Add columns for this species
    exon_result[[paste0(sp, "_TE")]] <- topTEs$TE_type
    exon_result[[paste0(sp, "_Count")]] <- topTEs$Count
    exon_result[[paste0(sp, "_pbCAZymes")]] <- topTEs$pbCAZymes
  }
  
  results_list[[etype]] <- exon_result
}

# -----------------------------
# CONCATENATE AND EXPORT
# -----------------------------

final_table <- bind_rows(results_list)

print(final_table)

# out_path <- file.path(base_path, "summary_Top4_TEs_with_size_per_species.tsv")
# fwrite(final_table, out_path, sep = "\t", quote = FALSE, col.names = TRUE)

# message("Final table saved as: ", out_path)
