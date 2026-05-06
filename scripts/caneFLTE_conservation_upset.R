# =========================================================
#
#                             Description
#
# Conservation analysis of full-length sugarcane TEs
# in grass genomes using RepeatMasker .out files
#
# Objectives:
# 1) Read the full-length sugarcane TE FASTA library and obtain the length of each TE_ID
# 2) Read all RepeatMasker .out files generated using the reclassified full-length sugarcane TE sequence library from TEdistill
# 3) Filter hits with coverage >= 80% of the full-length TE sequence
# 4) Select the best hit for each TE_ID in each genome
# 5) Build a complete presence/absence matrix, including fully absent sequences
# 6) Generate a summary table and an UpSet plot
#
# =========================================================

# -----------------------------
# 0. Parameters
# -----------------------------

root_dir <- "/home/mariacamiladiazrodriguez/ResultsTEs-MaCamiDR/RStudio/caneLib_in_grasses"
fasta_file <- file.path(root_dir, "distilledTE.flTE.renamed.fa")
coverage_threshold <- 0.80

best_hits_output_file <- "best_hits_caneLib_in_grasses_80pct.tsv"
presence_output_file <- "canesLib_ids_in_grasses_80pct_fullLength.tsv"
plot_output_file <- "UpSet_caneLib_in_grasses_80pct_fullLength.png"

# -----------------------------
# 1. Helper functions
# -----------------------------

# Remove parentheses and convert values to numeric
parse_rm_num <- function(x) {
  x %>%
    str_replace_all("[()]", "") %>%
    na_if("") %>%
    as.numeric()
}

# Read FASTA file and return TE_ID + sequence length
read_fasta_lengths <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  header_idx <- grep("^>", lines)

  if (length(header_idx) == 0) {
    stop("No FASTA headers found in: ", fasta_path)
  }

  end_idx <- c(header_idx[-1] - 1, length(lines))
  headers <- lines[header_idx]

  te_ids <- headers %>%
    str_remove("^>") %>%
    str_extract("^[^[:space:]]+") %>%
    str_remove("#.*$")

  te_lengths <- purrr::map2_int(header_idx, end_idx, function(i, j) {
    if (j <= i) return(0L)
    seq_lines <- lines[(i + 1):j]
    seq_lines <- str_replace_all(seq_lines, "\\s+", "")
    sum(nchar(seq_lines))
  })

  tibble(
    TE_ID = te_ids,
    TE_length = te_lengths
  ) %>%
    distinct(TE_ID, .keep_all = TRUE)
}

# Read a RepeatMasker .out file and return a processed table
read_repeatmasker_out <- function(path) {
  genome <- basename(path) %>%
    str_remove("_genomic\\.fna\\.out$")

  data <- read.table(
    path,
    skip = 3,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )

  if (ncol(data) < 15) {
    stop(".out file with fewer than 15 columns: ", path)
  }

  names(data) <- paste0("X", seq_len(ncol(data)))

  col_x15 <- if ("X15" %in% names(data)) data$X15 else rep(NA_character_, nrow(data))
  col_x16 <- if ("X16" %in% names(data)) data$X16 else rep(NA_character_, nrow(data))

  tibble(
    score = as.numeric(data$X1),
    perc_div = as.numeric(data$X2),
    perc_del = as.numeric(data$X3),
    perc_ins = as.numeric(data$X4),
    query_seq = data$X5,
    q_begin = as.numeric(data$X6),
    q_end = as.numeric(data$X7),
    q_left = data$X8,
    strand = data$X9,
    TE_ID = data$X10,
    class_family = data$X11,
    r_col12 = data$X12,
    r_col13 = data$X13,
    r_col14 = data$X14,
    rm_internal_id = col_x15,
    overlap_flag = col_x16,
    Genome = genome
  ) %>%
    mutate(
      query_aln_len = q_end - q_begin + 1,
      r12_num = parse_rm_num(r_col12),
      r13_num = parse_rm_num(r_col13),
      r14_num = parse_rm_num(r_col14),
      repeat_aln_len = case_when(
        str_detect(r_col12, "^\\(") ~ abs(r13_num - r14_num) + 1,
        TRUE ~ abs(r12_num - r13_num) + 1
      )
    ) %>%
    filter(
      !str_detect(
        class_family,
        regex("Simple_repeat|Low_complexity|Satellite|Unknown", ignore_case = TRUE)
      )
    )
}

# -----------------------------
# 2. Read the full-length library
# -----------------------------

te_library <- read_fasta_lengths(fasta_file)

if (nrow(te_library) == 0) {
  stop("The FASTA library was read, but no TE_ID was recovered.")
}

# -----------------------------
# 3. Locate .out files
# -----------------------------

out_files <- list.files(
  path = root_dir,
  pattern = "\\.out$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(out_files) == 0) {
  stop("No .out file was found in: ", root_dir)
}

genomes <- basename(out_files) %>%
  str_remove("_genomic\\.fna\\.out$") %>%
  unique()

if (length(genomes) == 0) {
  stop("Genome names could not be recovered from the .out files.")
}

# -----------------------------
# 4. Read and consolidate RepeatMasker hits
# -----------------------------

all_hits <- purrr::map_df(out_files, read_repeatmasker_out)

if (nrow(all_hits) == 0) {
  stop("After reading the .out files, no hit remained.")
}

# -----------------------------
# 5. Calculate coverage relative to the full-length TE sequence
# -----------------------------

all_hits <- all_hits %>%
  left_join(te_library, by = "TE_ID") %>%
  mutate(
    prop_full_length = repeat_aln_len / TE_length,
    perc_full_length = prop_full_length * 100
  )

te_without_length <- all_hits %>%
  filter(is.na(TE_length)) %>%
  distinct(TE_ID)

if (nrow(te_without_length) > 0) {
  warning(
    "There are TE_IDs in the .out files that were not found in the FASTA file. Examples: ",
    paste(head(te_without_length$TE_ID, 10), collapse = ", ")
  )
}

all_hits <- all_hits %>%
  filter(!is.na(TE_length), TE_length > 0)

# -----------------------------
# 6. Filter by minimum full-length coverage
# -----------------------------

approved_hits <- all_hits %>%
  filter(prop_full_length >= coverage_threshold)

if (nrow(approved_hits) == 0) {
  stop("No hit passed the coverage filter >= ", coverage_threshold * 100, "%.")
}

# -----------------------------
# 7. Select the best hit per TE_ID per genome
# -----------------------------

best_hit_per_genome <- approved_hits %>%
  arrange(
    Genome,
    TE_ID,
    desc(prop_full_length),
    perc_div,
    desc(score),
    desc(query_aln_len)
  ) %>%
  group_by(Genome, TE_ID) %>%
  slice(1) %>%
  ungroup()

write_tsv(best_hit_per_genome, best_hits_output_file)

# -----------------------------
# 8. Complete presence/absence matrix
# -----------------------------

presence_matrix_long <- expand_grid(
  TE_ID = te_library$TE_ID,
  Genome = genomes
) %>%
  left_join(
    best_hit_per_genome %>%
      transmute(TE_ID, Genome, Presence = 1L),
    by = c("TE_ID", "Genome")
  ) %>%
  mutate(Presence = replace_na(Presence, 0L))

presence_matrix <- presence_matrix_long %>%
  pivot_wider(
    names_from = Genome,
    values_from = Presence,
    values_fill = 0
  ) %>%
  left_join(te_library, by = "TE_ID") %>%
  relocate(TE_length, .after = TE_ID)

write_tsv(presence_matrix, presence_output_file)

# -----------------------------
# 9. Metadata by TE_ID for the UpSet plot
# -----------------------------

te_metadata <- best_hit_per_genome %>%
  group_by(TE_ID) %>%
  summarise(
    n_genomes_detected = n_distinct(Genome),
    mean_divergence_best_hits = mean(perc_div, na.rm = TRUE),
    min_divergence_best_hits = min(perc_div, na.rm = TRUE),
    mean_coverage_best_hits = mean(prop_full_length, na.rm = TRUE),
    max_coverage_best_hits = max(prop_full_length, na.rm = TRUE),
    .groups = "drop"
  )

upset_data <- presence_matrix %>%
  left_join(te_metadata, by = "TE_ID")

genome_columns <- genomes

# -----------------------------
# 10. Useful numerical summaries
# -----------------------------

genome_summary <- presence_matrix %>%
  summarise(across(all_of(genome_columns), sum)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Genome",
    values_to = "n_TEs_present"
  )

print(genome_summary)

total_library_te_ids <- nrow(te_library)

te_ids_detected_in_any_genome <- presence_matrix %>%
  filter(if_any(all_of(genome_columns), ~ .x == 1)) %>%
  nrow()

te_ids_absent_in_all_genomes <- presence_matrix %>%
  filter(if_all(all_of(genome_columns), ~ .x == 0)) %>%
  nrow()

cat("\n")
cat("Total TE_IDs in the library:", total_library_te_ids, "\n")
cat("TE_IDs with hit >= 80% in at least one genome:", te_ids_detected_in_any_genome, "\n")
cat("TE_IDs absent in all genomes:", te_ids_absent_in_all_genomes, "\n")
cat("\n")

# -----------------------------
# 11. UpSet plot
# -----------------------------

upset_plot <- upset(
  upset_data,
  genome_columns,
  name = "TE sharing",
  width_ratio = 0.20,
  annotations = list(
    "Mean divergence of best hits (%)" =
      ggplot(mapping = aes(y = mean_divergence_best_hits)) +
      geom_boxplot(fill = "lightblue", color = "darkblue", alpha = 0.7, na.rm = TRUE) +
      theme_minimal() +
      labs(y = "Divergence (%)", x = NULL)
  ),
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE)
  )
)

print(upset_plot)

# ggsave(filename = plot_output_file, plot = upset_plot, width = 11, height = 8.5, dpi = 300)
