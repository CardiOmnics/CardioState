# =============================================================================
# 9_revision_checks.R
# Checks reported in the manuscript text and Supplementary Table 5.
# Standalone: reads the newest *_good_paced cell tables written by the
# pipeline (well QC already applied), repeats the k-means clustering exactly
# as in 5_plot_and_stats.R, and writes one workbook:
#   3_output/Tables/revision_checks_<date>.xlsx
#
#   9   Nucleus area vs BODIPY (contraction session)
#   10  Cell area, MYH7 R403L vs MYH7 ISO
#   11  Cell area and N/C ratio, well-average model vs mixed model
#   12  Condition effect with and without k-means cluster in the mixed model
#   13  Oleate enrichment in contraction Cluster 3 (KW + Dunn, per line)
#   14  Nucleus area (S2) condition effect
#   15  Cluster composition per line and condition, change vs Control
#   16  Well and cell counts
# Run from 2_analysis/:  Rscript 9_revision_checks.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(lme4); library(lmerTest); library(emmeans)
  library(openxlsx); library(car)
})
source("../0_setup/colors.R")
source("../0_setup/helpers.R")
source("../0_setup/cluster_helpers.R")
emm_options(pbkrtest.limit = 50000)   # as in 5_plot_and_stats.R (Kenward-Roger df)

table_dir <- "../3_output/Tables"
newest <- function(stem) {
  f <- list.files(table_dir, pattern = paste0("^", stem, "_[0-9-]+\\.csv$"), full.names = TRUE)
  f[which.max(file.mtime(f))]
}
prep <- function(df) {
  relabel_condition(df) %>%
    # tables may carry the plate code ("1.6.24") or the table label
    # ("MYH7 ISO (1.6.24)"); keep only the code
    mutate(Cell_Line = sub("^.*\\((.*)\\)$", "\\1", as.character(Cell_Line)),
           Condition = factor(Condition, levels = group_levels),
           Cell_Line = factor(Cell_Line, levels = c("1.6", "1.6.24")))
}
ca  <- prep(read_csv(newest("calcium_good_paced"),     show_col_types = FALSE))
con <- prep(read_csv(newest("contraction_good_paced"), show_col_types = FALSE))
# Per-well QC (Supplementary Table 1, written by run_pipeline.R): all imaged wells
s1f <- list.files(table_dir, pattern = "^supplementary_table_1_[0-9-]+\\.xlsx$", full.names = TRUE)
s1x <- read.xlsx(s1f[which.max(file.mtime(s1f))], colNames = FALSE)
s1h <- which(s1x[, 1] == "modality")
wells_qc <- s1x[(s1h + 1):nrow(s1x), ]; names(wells_qc) <- unlist(s1x[s1h, ])
wells_qc$n_total <- as.numeric(wells_qc$n_total)
line_lab <- c("1.6" = "MYH7 R403L", "1.6.24" = "MYH7 ISO")
out <- list()
fmt_p <- function(p) format.pval(p, digits = 2, eps = 1e-16)

# ---- Clusters (same procedure as run_subpop_analysis) -----------------------
cluster_cells <- function(data, cols, label_rule, k = 3L) {
  cd <- data %>% select(cell_uid, well_uid, Condition, Cell_Line, all_of(cols)) %>% drop_na()
  m  <- cd %>% select(all_of(cols)) %>% scale()
  set.seed(42)
  km <- kmeans(m, centers = k, nstart = 50, iter.max = 200)
  km <- relabel_kmeans(km, m, label_rule)
  cd$cluster <- factor(km$cluster)
  cd
}
ca_cl <- cluster_cells(ca, c("ca_amp", "ca_CTD50_ms", "ca_CTD90_ms", "ca_rise_ms",
                             "ca_decay_ms", "Cell.Area", "Nucleus.Area", "NC_ratio"),
                       list(list(label = 3L, var = "ca_CTD50_ms", pick = "max"),
                            list(label = 2L, var = "ca_CTD50_ms", pick = "min")))
con_cl <- cluster_cells(con, c("contr_amp", "contr_rise_ms", "contr_decay_ms",
                               "Nucleus.Area", "Bodipy_intensity"),
                        list(list(label = 3L, var = "Bodipy_intensity", pick = "max"),
                             list(label = 1L, var = "contr_amp",        pick = "max")))
cat("Calcium cluster sizes:",     table(ca_cl$cluster),  "\n")
cat("Contraction cluster sizes:", table(con_cl$cluster), "\n")

# =============================================================================
# 9. Nucleus area vs BODIPY, contraction session
# =============================================================================
sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
nb <- con %>% filter(is.finite(Nucleus.Area), is.finite(Bodipy_intensity))
q9 <- bind_rows(
  sp(nb$Nucleus.Area, nb$Bodipy_intensity) %>% mutate(level = "Single cell", group = "All"),
  nb %>% group_by(Cell_Line) %>% group_modify(~ sp(.x$Nucleus.Area, .x$Bodipy_intensity)) %>%
    ungroup() %>% mutate(level = "Single cell", group = line_lab[as.character(Cell_Line)]) %>% select(-Cell_Line),
  nb %>% group_by(Cell_Line, Condition) %>% group_modify(~ sp(.x$Nucleus.Area, .x$Bodipy_intensity)) %>%
    ungroup() %>% mutate(level = "Single cell",
                         group = paste(line_lab[as.character(Cell_Line)], Condition)) %>%
    select(-Cell_Line, -Condition),
  # within-well: rank correlation after centring both on the well median
  nb %>% group_by(well_uid) %>%
    mutate(nuc_c = rank(Nucleus.Area) / n(), bod_c = rank(Bodipy_intensity) / n()) %>%
    ungroup() %>% { sp(.$nuc_c, .$bod_c) } %>% mutate(level = "Within well (ranks per well)", group = "All"),
  nb %>% group_by(well_uid) %>%
    summarise(nuc = mean(Nucleus.Area), bod = mean(Bodipy_intensity), .groups = "drop") %>%
    { sp(.$nuc, .$bod) } %>% mutate(level = "Well means", group = "All wells")
) %>% select(level, group, n, rho, p)
# mixed model: log10 BODIPY on log10 nucleus area, design terms + well
fit9 <- lmer(log10(Bodipy_intensity) ~ log10(Nucleus.Area) + Condition * Cell_Line + (1 | well_uid),
             data = nb)
co9 <- summary(fit9)$coefficients["log10(Nucleus.Area)", ]
q9 <- bind_rows(q9, tibble(level = "Mixed model (log-log slope, design + well)",
                           group = sprintf("slope = %.3f (SE %.3f), t(%.0f) = %.1f",
                                           co9["Estimate"], co9["Std. Error"], co9["df"], co9["t value"]),
                           n = nrow(nb), rho = NA, p = co9["Pr(>|t|)"]))
out$Q9_Nucleus_vs_BODIPY <- q9

# =============================================================================
# 10. Cell area, R403L vs ISO
# =============================================================================
fit_area <- lmer(log10(Cell.Area) ~ Condition * Cell_Line + (1 | well_uid),
                 data = ca %>% filter(Cell.Area > 0))
aov_area <- anova(fit_area)
em_line <- emmeans(fit_area, ~ Cell_Line)        # averaged over conditions
ratio_all <- as.data.frame(confint(contrast(em_line, "revpairwise"))) %>%
  mutate(across(c(estimate, lower.CL, upper.CL), ~ 10^.x))
p_all <- as.data.frame(contrast(em_line, "revpairwise"))$p.value
em_lc <- emmeans(fit_area, ~ Cell_Line | Condition)
by_cond <- as.data.frame(contrast(em_lc, "revpairwise")) %>%
  left_join(as.data.frame(confint(contrast(em_lc, "revpairwise"))) %>%
              select(Condition, lower.CL, upper.CL), by = "Condition") %>%
  mutate(ratio_ISO_vs_R403L = 10^estimate, lo = 10^lower.CL, hi = 10^upper.CL)
gm <- as.data.frame(emmeans(fit_area, ~ Cell_Line | Condition)) %>%
  mutate(geo_mean_um2 = 10^emmean) %>% select(Condition, Cell_Line, geo_mean_um2) %>%
  pivot_wider(names_from = Cell_Line, values_from = geo_mean_um2, names_prefix = "geo_mean_")
raw <- ca %>% group_by(Condition, Cell_Line) %>%
  summarise(n_cells = n(), n_wells = n_distinct(well_uid),
            mean = mean(Cell.Area), median = median(Cell.Area), .groups = "drop")
raw_line <- ca %>% group_by(Cell_Line) %>%
  summarise(n_cells = n(), n_wells = n_distinct(well_uid),
            mean = mean(Cell.Area), median = median(Cell.Area), .groups = "drop")
# well-average version of the line comparison
wa_area <- ca %>% group_by(well_uid, Condition, Cell_Line) %>%
  summarise(y = log10(mean(Cell.Area)), .groups = "drop")
old <- options(contrasts = c("contr.sum", "contr.poly"))
wa_aov <- car::Anova(lm(y ~ Condition * Cell_Line, data = wa_area), type = 3)
options(old)
q10 <- bind_rows(
  tibble(item = "Mixed model, Cell_Line main effect (log10 cell area)",
         value = sprintf("F(%d, %.1f) = %.2f, P = %s", aov_area["Cell_Line", "NumDF"],
                         aov_area["Cell_Line", "DenDF"], aov_area["Cell_Line", "F value"],
                         fmt_p(aov_area["Cell_Line", "Pr(>F)"]))),
  tibble(item = "Mixed model, Condition x Cell_Line",
         value = sprintf("F(%d, %.1f) = %.2f, P = %s", aov_area["Condition:Cell_Line", "NumDF"],
                         aov_area["Condition:Cell_Line", "DenDF"], aov_area["Condition:Cell_Line", "F value"],
                         fmt_p(aov_area["Condition:Cell_Line", "Pr(>F)"]))),
  tibble(item = "ISO / R403L geometric-mean ratio, averaged over conditions (95% CI)",
         value = sprintf("%.2f (%.2f to %.2f), P = %s", ratio_all$estimate, ratio_all$lower.CL,
                         ratio_all$upper.CL, fmt_p(p_all))),
  tibble(item = "Well-average two-way ANOVA, Cell_Line (log10 of well mean)",
         value = sprintf("F(1, %d) = %.2f, P = %s", wa_aov["Residuals", "Df"],
                         wa_aov["Cell_Line", "F value"], fmt_p(wa_aov["Cell_Line", "Pr(>F)"])))
)
out$Q10_CellArea_Line <- q10
out$Q10_CellArea_ByCondition <- by_cond %>%
  select(Condition, ratio_ISO_vs_R403L, lo, hi, t.ratio, df, p.value) %>%
  left_join(gm, by = "Condition")
out$Q10_CellArea_Raw <- bind_rows(raw, raw_line %>% mutate(Condition = "All")) %>%
  mutate(Cell_Line = line_lab[as.character(Cell_Line)])

# =============================================================================
# 11 + 14. Condition effect, well-average LM vs mixed model
# =============================================================================
compare_models <- function(df, col, log, label) {
  d <- df %>% select(well_uid, Condition, Cell_Line, y = all_of(col)) %>% drop_na()
  if (log) d <- d %>% filter(y > 0) %>% mutate(y = log10(y))
  fit <- lmer(y ~ Condition * Cell_Line + (1 | well_uid), data = d)
  a <- anova(fit)
  wa <- df %>% group_by(well_uid, Condition, Cell_Line) %>%
    summarise(y = mean(.data[[col]], na.rm = TRUE), .groups = "drop") %>% drop_na()
  if (log) wa <- wa %>% mutate(y = log10(y))
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  w <- car::Anova(lm(y ~ Condition * Cell_Line, data = wa), type = 3)
  options(old)
  tibble(parameter = label,
         mixed_condition = sprintf("F(3, %.1f) = %.2f, P = %s", a["Condition", "DenDF"],
                                   a["Condition", "F value"], fmt_p(a["Condition", "Pr(>F)"])),
         wellavg_condition = sprintf("F(3, %d) = %.2f, P = %s", w["Residuals", "Df"],
                                     w["Condition", "F value"], fmt_p(w["Condition", "Pr(>F)"])),
         mixed_line = sprintf("F(1, %.1f) = %.2f, P = %s", a["Cell_Line", "DenDF"],
                              a["Cell_Line", "F value"], fmt_p(a["Cell_Line", "Pr(>F)"])),
         wellavg_line = sprintf("F(1, %d) = %.2f, P = %s", w["Residuals", "Df"],
                                w["Cell_Line", "F value"], fmt_p(w["Cell_Line", "Pr(>F)"])),
         p_mixed = a["Condition", "Pr(>F)"], p_wellavg = w["Condition", "Pr(>F)"])
}
out$Q11_Q14_Models <- bind_rows(
  compare_models(ca,  "Cell.Area",    TRUE,  "Cell area"),
  compare_models(ca,  "NC_ratio",     TRUE,  "N/C ratio"),
  compare_models(ca,  "Nucleus.Area", TRUE,  "Nucleus area (S1)"),
  compare_models(con, "Nucleus.Area", TRUE,  "Nucleus area (S2)")
)
# post hoc for nucleus area S2 (Tukey within line, as in Supp. Table 3)
fit_n2 <- lmer(log10(Nucleus.Area) ~ Condition * Cell_Line + (1 | well_uid),
               data = con %>% filter(Nucleus.Area > 0))
out$Q14_NucS2_Posthoc <- as.data.frame(emmeans(fit_n2, pairwise ~ Condition | Cell_Line,
                                               adjust = "tukey")$contrasts) %>%
  mutate(pct_change = 100 * (10^(-estimate) - 1),   # second group relative to first
         Cell_Line = line_lab[as.character(Cell_Line)])
out$Q14_NucS2_Means <- as.data.frame(emmeans(fit_n2, ~ Condition | Cell_Line)) %>%
  mutate(geo_mean_um2 = 10^emmean, Cell_Line = line_lab[as.character(Cell_Line)]) %>%
  select(Cell_Line, Condition, geo_mean_um2)

# =============================================================================
# 12. Condition effect with vs without cluster membership
# =============================================================================
cluster_adjust <- function(cl_data, raw, col, log, label) {
  d <- cl_data %>% select(cell_uid, cluster) %>%
    inner_join(raw %>% select(cell_uid, well_uid, Condition, Cell_Line, y = all_of(col)),
               by = "cell_uid") %>% drop_na()
  if (log) d <- d %>% filter(y > 0) %>% mutate(y = log10(y))
  f0 <- lmer(y ~ Condition * Cell_Line + (1 | well_uid), data = d)
  f1 <- lmer(y ~ Condition * Cell_Line + cluster + (1 | well_uid), data = d)
  a0 <- anova(f0); a1 <- anova(f1)
  tibble(parameter = label, n_cells = nrow(d),
         F_without = a0["Condition", "F value"], p_without = a0["Condition", "Pr(>F)"],
         F_with_cluster = a1["Condition", "F value"], p_with_cluster = a1["Condition", "Pr(>F)"],
         p_line_without = a0["Cell_Line", "Pr(>F)"], p_line_with = a1["Cell_Line", "Pr(>F)"],
         p_int_without = a0["Condition:Cell_Line", "Pr(>F)"],
         p_int_with = a1["Condition:Cell_Line", "Pr(>F)"])
}
# Every parameter tested here was also a clustering input, so adjusting for
# cluster is partly circular. Read as exploratory.
q12 <- bind_rows(
  cluster_adjust(ca_cl, ca, "ca_amp",       TRUE,  "Calcium amplitude"),
  cluster_adjust(ca_cl, ca, "ca_CTD50_ms",  TRUE,  "Calcium CTD50"),
  cluster_adjust(ca_cl, ca, "ca_CTD90_ms",  TRUE,  "Calcium CTD90"),
  cluster_adjust(ca_cl, ca, "ca_rise_ms",   TRUE,  "Calcium rise time"),
  cluster_adjust(ca_cl, ca, "ca_decay_ms",  TRUE,  "Calcium decay time"),
  cluster_adjust(ca_cl, ca, "Cell.Area",    TRUE,  "Cell area"),
  cluster_adjust(ca_cl, ca, "Nucleus.Area", TRUE,  "Nucleus area (S1)"),
  cluster_adjust(ca_cl, ca, "NC_ratio",     TRUE,  "N/C ratio"),
  cluster_adjust(con_cl, con, "contr_amp",      TRUE, "Contraction amplitude"),
  cluster_adjust(con_cl, con, "contr_rise_ms",  TRUE, "Contraction rise time"),
  cluster_adjust(con_cl, con, "contr_decay_ms", TRUE, "Contraction decay time"),
  cluster_adjust(con_cl, con, "Nucleus.Area",   TRUE, "Nucleus area (S2)"),
  cluster_adjust(con_cl, con, "Bodipy_intensity", TRUE, "Bodipy intensity")
) %>%
  mutate(sig_without = p_without < 0.05, sig_with = p_with_cluster < 0.05,
         change = case_when(!sig_without & sig_with ~ "becomes significant",
                            sig_without & !sig_with ~ "loses significance",
                            p_with_cluster < p_without ~ "smaller P, same call",
                            TRUE ~ "larger P, same call"))
out$Q12_ClusterAdjusted <- q12

# =============================================================================
# 13 + 15. Cluster composition, KW and Dunn per line
# =============================================================================
well_pct <- function(cl) {
  cl %>% count(well_uid, Condition, Cell_Line, cluster) %>%
    complete(nesting(well_uid, Condition, Cell_Line), cluster, fill = list(n = 0)) %>%
    group_by(well_uid) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
}
kw_dunn <- function(wp, session) {
  wp %>% group_by(Cell_Line, cluster) %>% group_modify(function(d, g) {
    kt <- kruskal.test(pct ~ Condition, data = d)
    invisible(capture.output(
      dn <- FSA::dunnTest(pct ~ Condition, data = droplevels(d), method = "holm")$res))
    ctrl <- dn %>% filter(grepl("Control", Comparison))
    tibble(n_wells = n_distinct(d$well_uid),
           KW_chi_sq = unname(kt$statistic), KW_p = kt$p.value,
           Dunn_vs_Control = paste(sprintf("%s: P.adj = %.3f", ctrl$Comparison, ctrl$P.adj),
                                   collapse = "; "))
  }) %>% ungroup() %>% mutate(session = session, Cell_Line = line_lab[as.character(Cell_Line)],
                              .before = 1)
}
wp_ca  <- well_pct(ca_cl); wp_con <- well_pct(con_cl)
out$Q13_Q15_KW_Dunn <- bind_rows(kw_dunn(wp_ca, "Calcium"), kw_dunn(wp_con, "Contraction"))

comp <- function(cl, session) {
  cell_level <- cl %>% count(Cell_Line, Condition, cluster) %>%
    group_by(Cell_Line, Condition) %>% mutate(pct_cells = 100 * n / sum(n)) %>% ungroup()
  well_level <- well_pct(cl) %>% group_by(Cell_Line, Condition, cluster) %>%
    summarise(median_pct_wells = median(pct), n_wells = n(), .groups = "drop")
  cell_level %>% left_join(well_level, by = c("Cell_Line", "Condition", "cluster")) %>%
    group_by(Cell_Line, cluster) %>%
    mutate(change_vs_control_pp = pct_cells - pct_cells[Condition == "Control"]) %>%
    ungroup() %>%
    mutate(session = session, Cell_Line = line_lab[as.character(Cell_Line)], .before = 1)
}
out$Q15_Composition <- bind_rows(comp(ca_cl, "Calcium"), comp(con_cl, "Contraction"))
out$Q13_Cluster3_PerWell <- wp_con %>% filter(cluster == "3") %>%
  mutate(Cell_Line = line_lab[as.character(Cell_Line)]) %>% arrange(Cell_Line, Condition, well_uid)

# Cell-level sensitivity for item 13: logistic GLMM, Cluster 3 membership,
# Oleate vs Control, per line, well as random intercept
q13_glmm <- map_dfr(c("1.6", "1.6.24"), function(l) {
  d <- con_cl %>% filter(Cell_Line == l) %>% mutate(c3 = as.integer(cluster == "3"))
  f <- glmer(c3 ~ Condition + (1 | well_uid), data = d, family = binomial,
             control = glmerControl(optimizer = "bobyqa"))
  f0 <- update(f, . ~ . - Condition)
  co <- summary(f)$coefficients
  tibble(Cell_Line = line_lab[l],
         LRT_condition_p = anova(f0, f)$`Pr(>Chisq)`[2],
         oleate_vs_control_logOR = co["ConditionOleate", "Estimate"],
         oleate_vs_control_p = co["ConditionOleate", "Pr(>|z|)"],
         singular = isSingular(f))
})
out$Q13_GLMM_sensitivity <- q13_glmm

# =============================================================================
# 16. Well and cell counts
# =============================================================================
wells_ca_img  <- wells_qc$well_uid[wells_qc$modality == "Calcium"]
wells_con_img <- wells_qc$well_uid[wells_qc$modality == "Contraction"]
seg_retained  <- function(m) sum(wells_qc$n_total[wells_qc$modality == m & wells_qc$qc_status == "passed"])
counts <- tribble(
  ~item, ~value,
  "Calcium wells imaged", length(wells_ca_img),
  "Contraction wells imaged", length(wells_con_img),
  "Calcium wells retained (well QC passed)", n_distinct(ca$well_uid),
  "Contraction wells retained", n_distinct(con$well_uid),
  "Retained contraction wells that are also retained calcium wells",
    length(intersect(unique(con$well_uid), unique(ca$well_uid))),
  "Imaged contraction wells that are also imaged calcium wells",
    length(intersect(wells_con_img, wells_ca_img)),
  "Calcium electrically responsive cells (retained wells)", nrow(ca),
  "Calcium segmented cells, retained wells", seg_retained("Calcium"),
  "Calcium complete cases for clustering", nrow(ca_cl),
  "Calcium cells with CTD90 missing", sum(is.na(ca$ca_CTD90_ms)),
  "Contraction electrically responsive cells (retained wells)", nrow(con),
  "Contraction segmented cells, retained wells", seg_retained("Contraction"),
  "Contraction complete cases for clustering", nrow(con_cl)
)
# cell_uid is <plate>_<well>_<index within session>: equal IDs in the two
# sessions are NOT the same cell, so cells cannot be matched across sessions.
out$Q16_Counts <- counts
out$Q16_Wells_per_group <- bind_rows(
  ca  %>% distinct(well_uid, Condition, Cell_Line) %>% count(Cell_Line, Condition, name = "wells") %>% mutate(session = "Calcium"),
  con %>% distinct(well_uid, Condition, Cell_Line) %>% count(Cell_Line, Condition, name = "wells") %>% mutate(session = "Contraction")
) %>% mutate(Cell_Line = line_lab[as.character(Cell_Line)]) %>%
  pivot_wider(names_from = session, values_from = wells)

# ---- Write ---------------------------------------------------------------------
wb <- createWorkbook()
for (nm in names(out)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, out[[nm]])
  setColWidths(wb, nm, cols = seq_len(ncol(out[[nm]])), widths = "auto")
}
f_out <- file.path(table_dir, paste0("revision_checks_", format(Sys.Date(), "%Y-%m-%d"), ".xlsx"))
saveWorkbook(wb, f_out, overwrite = TRUE)
options(width = 200, pillar.sigfig = 4)
for (nm in names(out)) { cat("\n=====", nm, "=====\n"); print(as.data.frame(out[[nm]]), row.names = FALSE) }
cat("\nSaved:", f_out, "\n")
