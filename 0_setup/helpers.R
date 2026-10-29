library(dplyr)

prepare_groups <- function(df, group_col = "Group") {
  df %>%
    mutate(
      !!group_col := factor(
        .data[[group_col]],
        levels = group_levels
      )
    )
}

# -----------------------------------------------------------------------------
# Cell-line display labels
# -----------------------------------------------------------------------------
# The underlying Cell_Line values stay as the plate codes "1.6" / "1.6.24"
# (stable keys for joins and stats). These helpers provide the rendered
# strings used in figures and tables.

# Figures: italic MYH7 via plotmath expressions. Use with a facet labeller,
# e.g. facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller).
cell_line_fig_labels <- c(
  "1.6"    = "italic('MYH7')~'R403L'",
  "1.6.24" = "italic('MYH7')~'ISO'"
)
cell_line_fig_labeller <-
  ggplot2::as_labeller(cell_line_fig_labels, ggplot2::label_parsed)

# Tables: plain string with plate code retained in parentheses.
cell_line_table_labels <- c(
  "1.6"    = "MYH7 R403L (1.6)",
  "1.6.24" = "MYH7 ISO (1.6.24)"
)

# Rewrite the Cell_Line column of a data.frame to the table-display strings.
# Call right before writing a CSV / xlsx sheet.
relabel_cell_line_table <- function(df, col = "Cell_Line") {
  if (col %in% names(df)) {
    vals   <- as.character(df[[col]])
    mapped <- cell_line_table_labels[vals]
    mapped[is.na(mapped)] <- vals[is.na(mapped)]   # keep unknown codes as-is
    df[[col]] <- unname(mapped)
  }
  df
}

# -----------------------------------------------------------------------------
# Data loading functions
# -----------------------------------------------------------------------------

# Loading runs via an explicit manifest — one row per well, pointing at that
# well's Primary_DataTable.csv.
build_cyteseer_manifest <- function(cyteseer_dirs) {
  purrr::pmap_dfr(
    list(cyteseer_dirs$root_dir, cyteseer_dirs$plate_id),
    function(root_dir, plate_id) {
      if (!dir.exists(root_dir)) {
        warning("Directory does not exist: ", root_dir)
        return(NULL)
      }
      all_csv <- list.files(path = root_dir, pattern = "\\.csv$",
                            recursive = TRUE, full.names = TRUE)
      files <- all_csv[grepl("Primary_DataTable\\.csv$", basename(all_csv),
                             ignore.case = TRUE)]
      if (length(files) == 0) {
        warning("No Primary_DataTable.csv files found in: ", root_dir)
        return(NULL)
      }
      tibble::tibble(Plate = plate_id,
                     Well  = basename(dirname(files)),
                     file  = files)
    }
  )
}

read_cyteseer_manifest <- function(manifest, skip = 71) {
  purrr::pmap_dfr(
    list(manifest$file, manifest$Plate, manifest$Well),
    function(f, plate, well) {
      df <- read.csv(f, skip = skip, header = TRUE, fill = TRUE,
                     check.names = TRUE)
      df$Plate       <- plate
      df$Well        <- well
      df$source_file <- f
      df
    }
  )
}

# Read one CyteSeer results folder. Kept for standalone / interactive use;
# the pipeline goes through load_and_merge_cyteseer().
read_cyteseer_singlecell <- function(root_dir, plate_id, skip = 71) {
  manifest <- build_cyteseer_manifest(
    tibble::tibble(root_dir = root_dir, plate_id = plate_id)
  )
  if (nrow(manifest) == 0) return(NULL)
  read_cyteseer_manifest(manifest, skip = skip)
}

load_plate_maps <- function() {
  plate_map_16 <- read.delim("../1_data/Plate_Map_1.6.txt",
                             stringsAsFactors = FALSE) %>%
    mutate(Well      = as.character(Well),
           Condition = as.character(Condition),
           Cell_Line = as.character(Cell_Line),
           Plate     = "1.6")

  plate_map_1624 <- read.delim("../1_data/Plate_Map_1.6.24.txt",
                               stringsAsFactors = FALSE) %>%
    mutate(Well      = as.character(Well),
           Condition = as.character(Condition),
           Cell_Line = as.character(Cell_Line),
           Plate     = "1.6.24")

  bind_rows(
    plate_map_16   %>% select(Plate, Well, Condition, Cell_Line),
    plate_map_1624 %>% select(Plate, Well, Condition, Cell_Line)
  )
}

load_and_merge_cyteseer <- function(cyteseer_dirs, plate_map_all, skip = 71) {
  manifest <- build_cyteseer_manifest(cyteseer_dirs)

  singlecell_raw <- read_cyteseer_manifest(manifest, skip = skip)

  cat("  Raw rows loaded:", nrow(singlecell_raw), "\n")

  singlecell <- singlecell_raw %>%
    left_join(plate_map_all, by = c("Plate", "Well")) %>%
    mutate(
      Cell_Line = Plate,
      well_uid  = paste(Plate, Well, sep = "_")
    ) %>%
    group_by(well_uid) %>%
    mutate(
      cell_id  = as.character(row_number()),
      cell_uid = paste(well_uid, cell_id, sep = "_")
    ) %>%
    ungroup()

  if (exists("group_levels", envir = .GlobalEnv)) {
    singlecell$Condition <- factor(singlecell$Condition,
                                   levels = get("group_levels", envir = .GlobalEnv))
  }
  singlecell$Cell_Line <- factor(singlecell$Cell_Line, levels = c("1.6", "1.6.24"))

  cat("  Cells after merge:", nrow(singlecell), "\n")
  singlecell
}

# A duration longer than the pacing interval cannot be measured within one
# beat: the threshold crossing then lies in a neighbouring beat or at the end
# of the recording (typically when the signal does not fall back to the
# threshold between beats). Such values are set to missing, like other
# durations that could not be measured.
mask_unmeasurable_durations <- function(df, cols, max_ms = 1000, label = "") {
  for (cl in intersect(cols, names(df))) {
    bad <- is.finite(df[[cl]]) & df[[cl]] > max_ms
    if (any(bad))
      cat(sprintf("  -> %s %s: %d values > %d ms set to NA (all segmented cells)\n",
                  label, cl, sum(bad), max_ms))
    df[[cl]][bad] <- NA_real_
  }
  df
}

# ---- Manuscript display names ----------------------------------------------------
# Parameter names as they appear in the manuscript, written as plotmath so that
# Ca2+ gets a real superscript. Keys are the parameter labels used in the
# pipeline (model_params). Use pm_labels() as a scale labels function, e.g.
# scale_y_discrete(labels = pm_labels).
param_display_pm <- c(
  "Calcium amplitude"      = "Ca^{2*'+'}~amplitude",
  "Calcium CTD50"          = "CTD[50]",
  "Calcium CTD90"          = "CTD[90]",
  "Calcium rise time"      = "Ca^{2*'+'}~rise~time",
  "Calcium decay time"     = "Ca^{2*'+'}~decay~time",
  "Cell area"              = "Cell~area",
  "Nucleus area (S1)"      = "Nucleus~area~'(S1)'",
  "N/C ratio"              = "'N/C ratio'",
  "Contraction amplitude"  = "Contraction~amplitude",
  "Contraction rise time"  = "Contraction~rise~time",
  "Contraction decay time" = "Contraction~decay~time",
  "Nucleus area (S2)"      = "Nucleus~area~'(S2)'",
  "Bodipy intensity"       = "BODIPY~intensity"
)
pm_labels <- function(x) {
  x <- as.character(x)
  pm <- ifelse(x %in% names(param_display_pm), param_display_pm[x],
               paste0("'", gsub("'", "\\\\'", x), "'"))
  parse(text = pm)
}

# P value as plotmath, in the manuscript's style: P = 0.0032, P = 1.3 x 10^-9
fmt_p_pm <- function(p) {
  if (p >= 0.001) return(sprintf("italic(P)=='%s'", formatC(p, digits = 2, format = "fg", flag = "#")))
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("italic(P)=='%.1f'%%*%%10^%d", m, e)
}

