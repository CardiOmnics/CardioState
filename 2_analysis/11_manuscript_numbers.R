# =============================================================================
# 11_manuscript_numbers.R
# Every number the manuscript text cites, recomputed from the newest pipeline
# output, in the order it appears in the text. Used to check the Google Doc.
# Run from 2_analysis/ after run_pipeline.R:  Rscript 11_manuscript_numbers.R
# Output: printed list + 3_output/Tables/manuscript_numbers_<date>.csv
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(openxlsx) })
source("../0_setup/colors.R")
tdir <- "../3_output/Tables"
newest <- function(pat) { f <- list.files(tdir, pattern = pat, full.names = TRUE); f[which.max(file.mtime(f))] }
sr  <- newest("^statistical_results_[0-9-]+\\.xlsx$")
s1  <- newest("^supplementary_table_1_[0-9-]+\\.xlsx$")
sh  <- function(x) read.xlsx(sr, x)
ca  <- read_csv(newest("^calcium_good_paced_[0-9-]+\\.csv$"), show_col_types = FALSE)
con <- read_csv(newest("^contraction_good_paced_[0-9-]+\\.csv$"), show_col_types = FALSE)
cat("Source:", basename(sr), "\n")

out <- list()
add <- function(section, item, value) out[[length(out) + 1]] <<- tibble(section, item, value = as.character(value))
f3  <- function(x, d = 2) formatC(x, digits = d, format = "fg", flag = "#")
pf  <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, digits = 2, format = "fg"))
ln  <- c("1.6" = "R403L", "1.6.24" = "ISO")

an <- sh("ANOVA"); ph <- sh("Posthoc"); pc <- sh("Power_Comparison")
Fstat <- function(par, term) { r <- an[an$parameter == par & an$term == term, ]
  sprintf("F(%d, %.1f) = %.1f, P = %s", r$NumDF, r$DenDF, r$F.value, pf(r$`Pr(>F)`)) }
php <- function(par, line, contrast) { r <- ph[ph$parameter == par & grepl(line, ph$Cell_Line) & ph$contrast == contrast, ]
  sprintf("P = %s (estimate %.3f)", pf(r$p.value), r$estimate) }

# ---- Counts (Results, first section) ------------------------------------------
x <- read.xlsx(s1, colNames = FALSE); hdr <- which(x[, 1] == "modality")
w <- x[(hdr + 1):nrow(x), ]; names(w) <- unlist(x[hdr, ]); w$n_total <- as.numeric(w$n_total)
w$n_good_paced <- as.numeric(w$n_good_paced); w$n_irregular <- as.numeric(w$n_irregular)
for (m in c("Calcium", "Contraction")) { d <- w[w$modality == m, ]; p <- d$qc_status == "passed"
  add("Counts", paste(m, "wells imaged / retained"), paste(nrow(d), "/", sum(p)))
  add("Counts", paste(m, "% responsive, all imaged wells"), round(100 * sum(d$n_good_paced) / sum(d$n_total), 1))
  add("Counts", paste(m, "% irregular, retained wells"), round(100 * sum(d$n_irregular[p]) / sum(d$n_total[p]), 1)) }
add("Counts", "Calcium responsive cells / wells", paste(nrow(ca), "/", n_distinct(ca$well_uid)))
add("Counts", "Contraction responsive cells / wells", paste(nrow(con), "/", n_distinct(con$well_uid)))
wpg <- function(d) range(table(distinct(d, well_uid, Condition, Cell_Line)[, c("Condition", "Cell_Line")]))
add("Counts", "Wells per group, calcium", paste(wpg(ca), collapse = "-"))
add("Counts", "Wells per group, contraction", paste(wpg(con), collapse = "-"))
irr <- sh("Irregular_Summary")
ic <- irr[irr$session == "Contraction", ]
add("Counts", "Contraction % irregular per group (min, max)",
    paste(sprintf("%s %s %.0f%%", ic$Condition, ln[sub(".*\\((.*)\\)", "\\1", ic$Cell_Line)], ic$pct_irregular)[c(which.min(ic$pct_irregular), which.max(ic$pct_irregular))], collapse = "; "))
it <- sh("Irregular_Tests"); it <- it[it$modality == "Calcium", ]
add("Counts", "Irregular GLMM LRT calcium", paste(sprintf("%s chi2 = %.1f, P = %.3f", it$Cell_Line, it$LRT_chisq, it$LRT_p), collapse = "; "))

# ---- BODIPY ---------------------------------------------------------------------
sc <- sh("Summary_Contraction")
for (l in c("R403L", "ISO")) { r <- sc[grepl(l, sc$Cell_Line), ]
  add("BODIPY", paste("Mean fold oleate / control", l), round(r$Bodipy_intensity_mean[r$Condition == "Oleate"] / r$Bodipy_intensity_mean[r$Condition == "Control"], 2))
  add("BODIPY", paste("Mean fold oleate / palmitate", l), round(r$Bodipy_intensity_mean[r$Condition == "Oleate"] / r$Bodipy_intensity_mean[r$Condition == "Palmitate"], 2))
  for (cn in c("Control - Oleate", "Palmitate - Oleate", "Control - Palmitate", "Control - Combination"))
    add("BODIPY", paste(l, cn), php("Bodipy intensity", l, cn)) }

# ---- Ca2+ kinetics and morphology ----------------------------------------------------
for (p in c("Calcium amplitude", "Calcium CTD50", "Calcium CTD90", "Calcium rise time", "Calcium decay time",
            "Nucleus area (S1)", "Cell area", "N/C ratio"))
  add("Ca models", paste(p, "condition"), Fstat(p, "Condition"))
scal <- sh("Summary_Calcium")
r <- scal[grepl("R403L", scal$Cell_Line), ]
add("Ca models", "R403L mean decay (ms) Control/Palm/Ole/Comb", paste(round(r$ca_decay_ms_mean[match(group_levels, r$Condition)]), collapse = " / "))
add("Ca models", "R403L decay, combination vs control (%)", round(100 * (r$ca_decay_ms_mean[r$Condition == "Combination"] / r$ca_decay_ms_mean[r$Condition == "Control"] - 1)))
add("Ca models", "R403L decay Control - Combination", php("Calcium decay time", "R403L", "Control - Combination"))
add("Ca models", "R403L decay Control - Palmitate", php("Calcium decay time", "R403L", "Control - Palmitate"))
add("Ca models", "R403L rise, % change vs control (mean)", paste(round(100 * (r$ca_rise_ms_mean[match(c("Palmitate", "Oleate", "Combination"), r$Condition)] / r$ca_rise_ms_mean[r$Condition == "Control"] - 1)), collapse = " / "))
for (cn in c("Control - Palmitate", "Control - Oleate", "Control - Combination"))
  add("Ca models", paste("R403L rise", cn), php("Calcium rise time", "R403L", cn))
ri <- scal[grepl("ISO", scal$Cell_Line), ]
add("Ca models", "ISO decay, combination vs control (%)", round(100 * (ri$ca_decay_ms_mean[ri$Condition == "Combination"] / ri$ca_decay_ms_mean[ri$Condition == "Control"] - 1)))
add("Ca models", "ISO CTD50, combination vs control (%)", round(100 * (ri$ca_CTD50_ms_mean[ri$Condition == "Combination"] / ri$ca_CTD50_ms_mean[ri$Condition == "Control"] - 1)))
add("Ca models", "ISO decay Control - Combination", php("Calcium decay time", "ISO", "Control - Combination"))
add("Ca models", "ISO CTD50 Control - Combination", php("Calcium CTD50", "ISO", "Control - Combination"))
for (cn in c("Control - Palmitate", "Control - Oleate")) {
  add("Ca models", paste("ISO decay", cn), php("Calcium decay time", "ISO", cn))
  add("Ca models", paste("ISO CTD50", cn), php("Calcium CTD50", "ISO", cn)) }
add("Ca models", "Interaction rise", Fstat("Calcium rise time", "Condition:Cell_Line"))
add("Ca models", "Interaction decay", Fstat("Calcium decay time", "Condition:Cell_Line"))
for (l in c("R403L", "ISO")) { rr <- scal[grepl(l, scal$Cell_Line), ]
  add("Ca models", paste(l, "amplitude combination vs control (%)"), round(100 * (rr$ca_amp_mean[rr$Condition == "Combination"] / rr$ca_amp_mean[rr$Condition == "Control"] - 1)))
  add("Ca models", paste(l, "amplitude Control - Combination"), php("Calcium amplitude", l, "Control - Combination"))
  for (cn in c("Control - Palmitate", "Control - Oleate", "Control - Combination")) {
    e <- ph[ph$parameter == "Nucleus area (S1)" & grepl(l, ph$Cell_Line) & ph$contrast == cn, ]
    add("Ca models", paste(l, "nucleus S1", cn, "(% change, P)"), sprintf("%.0f%%, P = %s", 100 * (10^(-e$estimate) - 1), pf(e$p.value))) } }

# ---- Contraction ----------------------------------------------------------------------
for (p in c("Contraction rise time", "Contraction amplitude", "Contraction decay time", "Nucleus area (S2)"))
  add("Contraction models", paste(p, "condition"), Fstat(p, "Condition"))
for (cn in c("Control - Palmitate", "Control - Oleate", "Control - Combination"))
  add("Contraction models", paste("R403L rise", cn), php("Contraction rise time", "R403L", cn))
wc <- sh("Well_Cross_Cor")
add("Cross-modal", "n wells / rho range / smallest P", sprintf("%d / %.2f to %.2f / %s", wc$n_wells[1], min(wc$spearman_rho), max(wc$spearman_rho), pf(min(wc$p_value))))
add("Cross-modal", "Ca decay vs contraction decay", with(wc[wc$x_var == "ca_decay_ms" & wc$y_var == "contr_decay_ms", ], sprintf("rho = %.2f, P = %.3f", spearman_rho, p_value)))
add("Cross-modal", "Ca rise vs contraction rise", with(wc[wc$x_var == "ca_rise_ms" & wc$y_var == "contr_rise_ms", ], sprintf("rho = %.2f, P = %.3f", spearman_rho, p_value)))

# ---- Line differences ------------------------------------------------------------------
for (p in c("Calcium decay time", "Calcium CTD50", "Calcium CTD90", "Contraction amplitude", "Cell area", "N/C ratio", "Nucleus area (S1)"))
  add("Line effects", p, Fstat(p, "Cell_Line"))
rc <- sc[sc$Condition == "Control", ]
add("Line effects", "Contraction amplitude control means R403L / ISO", paste(signif(rc$contr_amp_mean[grepl("R403L", rc$Cell_Line)], 2), "/", signif(rc$contr_amp_mean[grepl("ISO", rc$Cell_Line)], 2)))

# ---- Variance decomposition ----------------------------------------------------------------
vd <- sh("Variance_Decomp")
names(vd) <- c("parameter", "line", "condition", "well", "within")[seq_len(ncol(vd))]
add("Variance", "Parameters where within-well is largest", sprintf("%d of %d", sum(vd$within >= pmax(vd$line, vd$condition, vd$well)), nrow(vd)))
for (i in seq_len(nrow(vd))) add("Variance", vd$parameter[i], sprintf("line %.1f, condition %.1f, well %.1f, within %.1f", vd$line[i], vd$condition[i], vd$well[i], vd$within[i]))

# ---- Clusters ------------------------------------------------------------------------------
kd <- sh("K_Diagnostics")
for (s in unique(kd$session)) add("Clusters", paste(s, "silhouette k = 2..6"), paste(sprintf("%.2f", kd$mean_sil[kd$session == s]), collapse = " / "))
for (s in c("Ca", "Contr")) {
  pr <- sh(paste0("Clust_", s, "_Profiles"))
  for (i in seq_len(nrow(pr))) add("Clusters", paste(s, "cluster", pr$cluster[i], "profile"),
    paste(paste(names(pr)[1:(ncol(pr) - 2)], signif(unlist(pr[i, 1:(ncol(pr) - 2)]), 3), sep = " = "), collapse = ", ") %>% paste0(", n = ", pr$n_cells[i], " (", round(100 * pr$n_cells[i] / sum(pr$n_cells), 1), "%)"))
  cp <- sh(paste0("Clust_", s, "_Composition"))
  for (l in c("R403L", "ISO")) for (cd in group_levels) { z <- cp[grepl(l, cp$Cell_Line) & cp$Condition == cd, ]
    add("Clusters", paste(s, l, cd, "% in C1/C2/C3"), paste(round(z$pct[match(1:3, z$cluster)]), collapse = " / ")) }
  kw <- sh(paste0("Clust_", s, "_WellKW"))
  for (i in seq_len(nrow(kw))) add("Clusters", paste(s, "KW", kw$Cell_Line[i], "cluster", kw$cluster[i]), sprintf("chi2 = %.1f, P = %.4f, n = %d wells", kw$KW_chi_sq[i], kw$p_value[i], kw$n_wells[i]))
}

# ---- Fig. 6 ------------------------------------------------------------------------------
add("Fig 6", "Significant in both", paste(pc$parameter[pc$sig_wellavg == "*" & pc$sig_mixed == "*"], collapse = ", "))
add("Fig 6", "Non-significant in both", paste(pc$parameter[pc$sig_wellavg == "ns" & pc$sig_mixed == "ns"], collapse = ", "))
add("Fig 6", "Disagreements", paste(pc$parameter[pc$sig_wellavg != pc$sig_mixed], collapse = ", "))
wd <- sh("WithinCluster_Direction")
for (p in c("Cell area", "N/C ratio")) for (k in 1:3) { z <- wd[wd$parameter == p & wd$session == "calcium" & wd$cluster == k & wd$Condition != "Control", ]
  add("Fig 6 within-cluster", paste(p, "cluster", k), paste0(sprintf("P(condition) = %s; ", pf(z$p_condition[1])),
      paste(sprintf("%s %+.0f%% (P = %s)", z$Condition, z$pct_vs_Control, pf(z$p_vs_Control)), collapse = "; "))) }
sig_wc <- wd %>% distinct(parameter, session, cluster, p_condition) %>% group_by(parameter) %>%
  summarise(sig_clusters = paste(cluster[p_condition < 0.05], collapse = ","), .groups = "drop")
add("Fig 6 within-cluster", "Clusters with P < 0.05 per parameter", paste(sprintf("%s [%s]", sig_wc$parameter, sig_wc$sig_clusters), collapse = "; "))

res <- bind_rows(out)
options(width = 250)
for (s in unique(res$section)) { cat("\n=====", s, "=====\n"); r <- res[res$section == s, ]; for (i in seq_len(nrow(r))) cat(sprintf("  %-55s %s\n", r$item[i], r$value[i])) }
write_csv(res, file.path(tdir, paste0("manuscript_numbers_", format(Sys.Date(), "%Y-%m-%d"), ".csv")))
