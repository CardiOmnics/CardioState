# =============================================================================
# 10_supplementary_tables.R
# Builds the supplementary tables in their manuscript order, one sheet each,
# from the newest pipeline output (run_pipeline.R) and 9_revision_checks.R.
#
#   Supplementary Table 1  Well-level quality control
#   Supplementary Table 2  Cluster-number diagnostics (k = 2-6)
#   Supplementary Table 3  Linear mixed-model ANOVA
#   Supplementary Table 4  Post hoc contrasts (Tukey, within cell line)
#   Supplementary Table 5  Cluster profiles, composition and well-level tests
#   Supplementary Table 6  Condition effects within clusters
#
# Run from 2_analysis/:  Rscript 10_supplementary_tables.R
# Output: 3_output/Tables/Supplementary_Tables_<date>.xlsx
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(openxlsx) })
tdir <- "../3_output/Tables"
newest <- function(pat) { f <- list.files(tdir, pattern = pat, full.names = TRUE); f[which.max(file.mtime(f))] }
sr <- newest("^statistical_results_[0-9-]+\\.xlsx$")
s1 <- newest("^supplementary_table_1_[0-9-]+\\.xlsx$")
rc <- newest("^revision_checks_[0-9-]+\\.xlsx$")
cat("Using:", basename(sr), basename(s1), basename(rc), "\n")
sh <- function(x) read.xlsx(sr, x)

line_lab <- function(x) sub(" \\(.*\\)$", "", x)          # "MYH7 ISO (1.6.24)" -> "MYH7 ISO"
param_order <- c("Calcium amplitude", "Calcium CTD50", "Calcium CTD90", "Calcium rise time",
                 "Calcium decay time", "Cell area", "Nucleus area (S1)", "N/C ratio",
                 "Contraction amplitude", "Contraction rise time", "Contraction decay time",
                 "Nucleus area (S2)", "Bodipy intensity")
session_of <- function(p) ifelse(p %in% param_order[1:8], "Calcium (S1)", "Contraction (S2)")
cond_levels <- c("Control", "Palmitate", "Oleate", "Combination")

tables <- list()

# ---- 1. Well QC ----------------------------------------------------------------------
x <- read.xlsx(s1, colNames = FALSE)
hdr <- which(x[, 1] == "modality")
t1 <- x[(hdr + 1):nrow(x), ]; names(t1) <- unlist(x[hdr, ])
t1 <- t1 %>% transmute(Session = modality, Well = well_uid, `Cell line` = line_lab(Cell_Line),
                       Condition, `Segmented cells` = as.integer(n_total),
                       `Electrically responsive` = as.integer(n_good_paced),
                       Irregular = as.integer(n_irregular),
                       `% responsive` = as.numeric(pct_good_paced), `% irregular` = as.numeric(pct_irregular),
                       `% CyteSeer responders` = suppressWarnings(as.numeric(pct_responder)),
                       `QC status` = qc_status, `Reason for exclusion` = qc_fail_reason)
tables[["Supplementary Table 1"]] <- list(
  title = "Supplementary Table 1. Well-level quality control.",
  legend = paste("For every imaged well, the number of segmented cells, the number and percentage of electrically",
                 "responsive and irregular cells, the percentage of cells flagged as responders by CyteSeer (calcium",
                 "session only), the QC status and the reason for exclusion. Wells were excluded if they contained fewer",
                 "than 10 electrically responsive cells or, in the calcium session, if fewer than 71.0% of cells were",
                 "CyteSeer responders (Otsu threshold). Wells that beat at twice the pacing frequency were excluded by hand.",
                 "Abbreviations: QC, quality control."),
  data = t1)

# ---- 2. k diagnostics -------------------------------------------------------------------
t2 <- sh("K_Diagnostics") %>% transmute(Session = session, k, `Within-cluster sum of squares` = round(wss),
                                         `Mean silhouette width` = round(mean_sil, 3),
                                         `Cells clustered` = n_cells, `Cells in silhouette subsample` = n_sil_used)
tables[["Supplementary Table 2"]] <- list(
  title = "Supplementary Table 2. Cluster-number diagnostics.",
  legend = paste("Within-cluster sum of squares and mean silhouette width for k-means solutions with k = 2 to 6, per",
                 "imaging session, on z-scored clustering variables. Silhouette widths were computed on a random",
                 "subsample of 5,000 cells. Both sessions were clustered with k = 3.",
                 "Abbreviations: k, number of clusters."),
  data = t2)

# ---- 3. ANOVA ----------------------------------------------------------------------------
pc <- sh("Power_Comparison")
t3 <- sh("ANOVA") %>%
  left_join(pc %>% select(parameter, n_wells, n_cells), by = "parameter") %>%
  transmute(Session = session_of(parameter), Parameter = parameter,
            Term = dplyr::recode(term, Cell_Line = "Cell line", `Condition:Cell_Line` = "Condition x cell line"),
            `Sum of squares` = signif(Sum.Sq, 4), `Mean square` = signif(Mean.Sq, 4),
            `Numerator df` = NumDF, `Denominator df` = round(DenDF, 1), F = round(F.value, 2),
            P = signif(`Pr(>F)`, 3), Scale = ifelse(log_transformed, "log10", "untransformed"),
            Cells = n_cells, Wells = n_wells) %>%
  arrange(match(Parameter, param_order))
tables[["Supplementary Table 3"]] <- list(
  title = "Supplementary Table 3. Linear mixed-model ANOVA.",
  legend = paste("Type III F tests (Satterthwaite denominator degrees of freedom) for condition, cell line and",
                 "their interaction in the model parameter ~ condition \u00d7 cell line + (1 | well), fitted on",
                 "log10-transformed values of each parameter. Durations longer than the 1,000 ms pacing interval",
                 "were treated as missing. Abbreviations: df, degrees of freedom."),
  data = t3)

# ---- 4. Post hoc -------------------------------------------------------------------------
t4 <- sh("Posthoc") %>%
  transmute(Session = session_of(parameter), Parameter = parameter, `Cell line` = line_lab(Cell_Line),
            Contrast = contrast, `Estimate (log10)` = signif(estimate, 3),
            `Ratio (first / second)` = signif(10^estimate, 3), SE = signif(SE, 3),
            df = round(df, 1), t = round(t.ratio, 2), `Adjusted P` = signif(p.value, 3)) %>%
  arrange(match(Parameter, param_order), `Cell line`)
tables[["Supplementary Table 4"]] <- list(
  title = "Supplementary Table 4. Post hoc pairwise comparisons.",
  legend = paste("Pairwise comparisons between conditions within each cell line from the models in Supplementary",
                 "Table 3 (estimated marginal means, Tukey adjustment, Kenward-Roger degrees of freedom). Estimates",
                 "are differences on the log10 scale, and the ratio gives the first condition relative to the second.",
                 "Abbreviations: SE, standard error; df, degrees of freedom."),
  data = t4)

# ---- 5. Clusters ---------------------------------------------------------------------------
prof <- bind_rows(
  sh("Clust_Ca_Profiles") %>% mutate(Session = "Calcium (S1)"),
  sh("Clust_Contr_Profiles") %>% mutate(Session = "Contraction (S2)")) %>%
  rename(Cluster = cluster, Cells = n_cells) %>%
  rename_with(~ dplyr::recode(.x, ca_amp = "Ca amplitude (a.u.)", ca_CTD50_ms = "CTD50 (ms)", ca_CTD90_ms = "CTD90 (ms)",
                       ca_rise_ms = "Ca rise time (ms)", ca_decay_ms = "Ca decay time (ms)",
                       Cell.Area = "Cell area (um2)", Nucleus.Area = "Nucleus area (um2)", NC_ratio = "N/C ratio",
                       contr_amp = "Contraction amplitude (a.u.)", contr_rise_ms = "Contraction rise time (ms)",
                       contr_decay_ms = "Contraction decay time (ms)", Bodipy_intensity = "BODIPY intensity (a.u.)")) %>%
  mutate(across(where(is.numeric) & !c(Cluster, Cells), ~ signif(.x, 3))) %>%
  relocate(Session, Cluster, Cells)
comp <- read.xlsx(rc, "Q15_Composition") %>%
  transmute(Session = paste0(session, ifelse(session == "Calcium", " (S1)", " (S2)")), `Cell line` = Cell_Line,
            Condition = factor(Condition, levels = cond_levels), Cluster = cluster, Cells = n,
            `% of cells` = round(pct_cells, 1), `Median % per well` = round(median_pct_wells, 1),
            Wells = n_wells, `Change vs control (percentage points)` = round(change_vs_control_pp, 1)) %>%
  # clusters without cells in a group are missing from the counts; show them as 0
  complete(nesting(Session, `Cell line`, Condition), Cluster,
           fill = list(Cells = 0L, `% of cells` = 0, `Median % per well` = 0)) %>%
  group_by(Session, `Cell line`, Condition) %>% fill(Wells, .direction = "downup") %>% ungroup() %>%
  group_by(Session, `Cell line`, Cluster) %>%
  mutate(`Change vs control (percentage points)` = round(`% of cells` - `% of cells`[Condition == "Control"], 1)) %>%
  ungroup() %>%
  arrange(Session, `Cell line`, Condition, Cluster)
kw <- read.xlsx(rc, "Q13_Q15_KW_Dunn") %>%
  transmute(Session = paste0(session, ifelse(session == "Calcium", " (S1)", " (S2)")), `Cell line` = Cell_Line,
            Cluster = cluster, Wells = n_wells, `Kruskal-Wallis chi2` = round(KW_chi_sq, 2), `KW P` = signif(KW_p, 3),
            `Dunn vs control (Holm-adjusted P)` = Dunn_vs_Control)
tables[["Supplementary Table 5"]] <- list(
  title = "Supplementary Table 5. Functional subpopulations: profiles, composition and well-level tests.",
  legend = paste("(A) Cluster centres in original units and cluster sizes. (B) Number and percentage of cells per",
                 "cluster for each condition and cell line, the median percentage across wells and the change from",
                 "control. (C) Kruskal-Wallis tests of the per-well cluster percentage across conditions, per cell line",
                 "and cluster, with Dunn post hoc tests against control (Holm adjustment). Calcium clusters are numbered",
                 "by kinetics (1 intermediate, 2 fastest, 3 slowest), contraction clusters by amplitude and BODIPY",
                 "(1 highest amplitude, 3 highest BODIPY). Abbreviations: KW, Kruskal-Wallis."),
  parts = list(A = prof, B = comp, C = kw))

# ---- 6. Within-cluster effects -----------------------------------------------------------------
t6 <- sh("WithinCluster_Direction") %>%
  transmute(Session = ifelse(session == "calcium", "Calcium (S1)", "Contraction (S2)"), Parameter = parameter,
            Cluster = cluster, Cells = n_cells, Wells = n_wells, `P (condition)` = signif(p_condition, 3),
            Method = dplyr::recode(method, lmer = "mixed model", `KW fallback` = "Kruskal-Wallis on well means"),
            Condition = factor(Condition, levels = cond_levels), `Mean (original units)` = signif(mean, 3),
            `% vs control` = round(pct_vs_Control, 1), `Dunnett P vs control` = signif(p_vs_Control, 3)) %>%
  arrange(match(Parameter, param_order), Cluster, Condition)
tables[["Supplementary Table 6"]] <- list(
  title = "Supplementary Table 6. Condition effects within clusters.",
  legend = paste("Condition effect tested separately within each cluster with parameter ~ condition + cell line +",
                 "(1 | well) on log10-transformed values (Kruskal-Wallis on well means where the mixed model could",
                 "not be estimated). Means are estimated marginal means back-transformed to original units, averaged",
                 "over cell line, with Dunnett-adjusted contrasts against control. Because the parameters were also",
                 "clustering variables, these tests are exploratory."),
  data = t6)

# ---- Write -------------------------------------------------------------------------------------
wb <- createWorkbook()
st_title <- createStyle(textDecoration = "bold", fontSize = 12)
st_leg   <- createStyle(wrapText = TRUE, valign = "top")
st_head  <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom", wrapText = TRUE)
st_part  <- createStyle(textDecoration = "bold", fontColour = "#305496")
for (nm in names(tables)) {
  t <- tables[[nm]]; addWorksheet(wb, nm)
  writeData(wb, nm, t$title, startRow = 1); addStyle(wb, nm, st_title, rows = 1, cols = 1)
  writeData(wb, nm, t$legend, startRow = 2); addStyle(wb, nm, st_leg, rows = 2, cols = 1)
  parts <- if (!is.null(t$parts)) t$parts else list(` ` = t$data)
  # CTD50 / CTD90 with Unicode subscript digits, in column names and cell text
  sub_ctd <- function(x) gsub("CTD90", "CTD\u2089\u2080", gsub("CTD50", "CTD\u2085\u2080", x))
  parts <- lapply(parts, function(d) { names(d) <- sub_ctd(names(d))
    dplyr::mutate(d, dplyr::across(where(is.character), sub_ctd)) })
  ncols <- max(sapply(parts, ncol)); mergeCells(wb, nm, cols = 1:ncols, rows = 2); setRowHeights(wb, nm, 2, 75)
  row <- 4
  for (pn in names(parts)) {
    if (pn != " ") { writeData(wb, nm, paste0("(", pn, ")"), startRow = row); addStyle(wb, nm, st_part, rows = row, cols = 1); row <- row + 1 }
    writeData(wb, nm, parts[[pn]], startRow = row, headerStyle = st_head)
    row <- row + nrow(parts[[pn]]) + 3
  }
  setColWidths(wb, nm, cols = 1:ncols, widths = 16)
}
out <- file.path(tdir, paste0("Supplementary_Tables_", format(Sys.Date(), "%Y-%m-%d"), ".xlsx"))
saveWorkbook(wb, out, overwrite = TRUE)
cat("Saved:", out, "\n")
for (nm in names(tables)) { t <- tables[[nm]]; n <- if (!is.null(t$parts)) paste(sapply(t$parts, nrow), collapse = "+") else nrow(t$data); cat(" ", nm, ":", n, "rows\n") }
