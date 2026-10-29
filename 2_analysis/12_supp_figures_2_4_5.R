# =============================================================================
# 12_supp_figures_2_4_5.R
# Publication versions (180 mm wide, text 5-7 pt) of three supplementary
# figures that 5_plot_and_stats.R draws as exploratory plots:
#
#   Supplementary Fig. 2  Irregular cells per well, per session and line
#   Supplementary Fig. 4  Single-cell Spearman correlations (A calcium, B contraction)
#   Supplementary Fig. 5  Well-level correlations between Ca2+ and contraction
#
# Same data and statistics as the exploratory versions: per-well values from
# Supplementary Table 1, complete-case Spearman correlations on the
# electrically responsive cells (values shown where P < 0.05), and Spearman
# tests on well means for the seven reported pairs.
#
# Run from 2_analysis/ after run_pipeline.R:  Rscript 12_supp_figures_2_4_5.R
# Output: 3_output/Figures/Publication/supp_figures_2_4_5/
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(openxlsx); library(ggcorrplot) })
source("../0_setup/colors.R")
source("../0_setup/theme.R")
source("../0_setup/helpers.R")

main_dir      <- "../3_output"
analysis_date <- format(Sys.Date(), "%Y-%m-%d")
out_dir       <- file.path(main_dir, "Figures", "Publication", "supp_figures_2_4_5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tdir   <- file.path(main_dir, "Tables")
newest <- function(pat) { f <- list.files(tdir, pattern = pat, full.names = TRUE); f[which.max(file.mtime(f))] }
code   <- function(x) sub("^.*\\((.*)\\)$", "\\1", as.character(x))   # "MYH7 ISO (1.6.24)" -> "1.6.24"
prep   <- function(d) relabel_condition(d) %>%
  mutate(Cell_Line = factor(code(Cell_Line), levels = c("1.6", "1.6.24")),
         Condition = factor(Condition, levels = group_levels))

ca  <- prep(read_csv(newest("^calcium_good_paced_[0-9-]+\\.csv$"), show_col_types = FALSE))
con <- prep(read_csv(newest("^contraction_good_paced_[0-9-]+\\.csv$"), show_col_types = FALSE))

line_names <- c("1.6" = expression(italic('MYH7')~R403L), "1.6.24" = expression(italic('MYH7')~ISO))
line_lab   <- as_labeller(c("1.6" = "italic('MYH7')~R403L", "1.6.24" = "italic('MYH7')~ISO"), label_parsed)

# Base theme: 7 pt titles, 6 pt axis text, 5.5 pt small print
pub <- theme_bw(base_size = 7) +
  theme(axis.text = element_text(size = 6, colour = "grey20"), axis.title = element_text(size = 7),
        strip.text = element_text(size = 7), strip.background = element_blank(),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.4, fill = NA),
        axis.ticks = element_line(linewidth = 0.3), legend.key.size = unit(3, "mm"),
        plot.margin = margin(1.5, 2, 1.5, 1.5, "mm"))

# ---- Supplementary Fig. 2: irregular cells per well ------------------------------------------
s1x <- read.xlsx(newest("^supplementary_table_1_[0-9-]+\\.xlsx$"), colNames = FALSE)
h   <- which(s1x[, 1] == "modality")
wq  <- s1x[(h + 1):nrow(s1x), ]; names(wq) <- unlist(s1x[h, ])
wq  <- wq %>% relabel_condition() %>%
  transmute(Session = factor(dplyr::recode(modality, Calcium = "Calcium session (S1)", Contraction = "Contraction session (S2)"),
                             levels = c("Calcium session (S1)", "Contraction session (S2)")),
            Cell_Line = factor(code(Cell_Line), levels = c("1.6", "1.6.24")),
            Condition = factor(Condition, levels = group_levels),
            pct = as.numeric(pct_irregular),
            Status = ifelse(qc_status == "passed", "Included", "Excluded (well QC)"))
set.seed(1)
wq$x <- as.numeric(wq$Condition) + runif(nrow(wq), -0.22, 0.22)
p2 <- ggplot(wq, aes(x, pct)) +
  geom_point(data = filter(wq, Status == "Included"), aes(fill = Condition),
             shape = 21, size = 1.6, colour = "grey25", stroke = 0.25) +
  geom_point(data = filter(wq, Status != "Included"), shape = 17, size = 1.6, colour = "#D62728") +
  facet_grid(Session ~ Cell_Line, labeller = labeller(Cell_Line = line_lab)) +
  scale_x_continuous(breaks = 1:4, labels = group_levels) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Irregular cells (% of segmented cells)") +
  pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, paste0("supp_figure2_", analysis_date, ".pdf")), p2,
       width = 180, height = 95, units = "mm", device = cairo_pdf)

# ---- Supplementary Fig. 4: single-cell Spearman correlations -----------------------------------
ca_vars  <- c(ca_amp = "Ca^{2*'+'}~amplitude", ca_CTD50_ms = "CTD[50]", ca_CTD90_ms = "CTD[90]",
              ca_rise_ms = "Ca^{2*'+'}~rise~time", ca_decay_ms = "Ca^{2*'+'}~decay~time",
              Cell.Area = "Cell~area", Nucleus.Area = "Nucleus~area", NC_ratio = "'N/C ratio'")
con_vars <- c(contr_amp = "Contraction~amplitude", contr_rise_ms = "Contraction~rise~time",
              contr_decay_ms = "Contraction~decay~time", Nucleus.Area = "Nucleus~area",
              Bodipy_intensity = "BODIPY~intensity")
cor_long <- function(d, vars) {
  m <- d %>% select(all_of(names(vars))) %>% drop_na()
  r <- cor(m, method = "spearman"); p <- cor_pmat(m, method = "spearman")
  idx <- which(lower.tri(r), arr.ind = TRUE)
  list(n = nrow(m), d = tibble(v1 = factor(names(vars)[idx[, 1]], levels = names(vars)),
                               v2 = factor(names(vars)[idx[, 2]], levels = names(vars)),
                               rho = r[idx], p = p[idx]))
}
heat <- function(cl, vars, title) {
  ggplot(cl$d, aes(v2, v1, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(p < 0.05, sprintf("%.2f", rho), "")), size = 5 / .pt) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-1, 1),
                         name = expression("Spearman"~rho)) +
    scale_x_discrete(labels = function(x) parse(text = vars[x])) +
    scale_y_discrete(labels = function(x) parse(text = vars[x])) +
    labs(x = NULL, y = NULL, subtitle = sprintf("%s, n = %s cells", title, format(cl$n, big.mark = ","))) +
    coord_fixed() + pub +
    theme(axis.text = element_text(size = 7, colour = "grey20"),
          axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
          panel.border = element_blank(), plot.subtitle = element_text(size = 7, face = "bold"),
          legend.position = "bottom", legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm"))
}
c_ca <- cor_long(ca, ca_vars); c_co <- cor_long(con, con_vars)
p4 <- cowplot::plot_grid(heat(c_ca, ca_vars, "Calcium session (S1)"),
                         heat(c_co, con_vars, "Contraction session (S2)"),
                         ncol = 2, rel_widths = c(1.15, 1), labels = c("a", "b"), label_size = 8)
ggsave(file.path(out_dir, paste0("supp_figure4_", analysis_date, ".pdf")), p4,
       width = 180, height = 95, units = "mm", device = cairo_pdf)

# ---- Supplementary Fig. 5: well-level Ca2+ vs contraction --------------------------------------
wm <- inner_join(
  ca %>% group_by(well_uid, Condition, Cell_Line) %>%
    summarise(across(c(ca_amp, ca_CTD50_ms, ca_CTD90_ms, ca_rise_ms, ca_decay_ms), ~ mean(.x, na.rm = TRUE)), .groups = "drop"),
  con %>% group_by(well_uid) %>%
    summarise(across(c(contr_amp, contr_rise_ms, contr_decay_ms), ~ mean(.x, na.rm = TRUE)), .groups = "drop"),
  by = "well_uid")
pairs <- tribble(~x, ~y,
  "ca_decay_ms", "contr_amp", "ca_amp", "contr_amp", "ca_CTD50_ms", "contr_amp", "ca_CTD90_ms", "contr_amp",
  "ca_rise_ms", "contr_rise_ms", "ca_decay_ms", "contr_decay_ms", "ca_amp", "contr_decay_ms")
# units quoted, so plotmath does not draw the brackets as scaled delimiters
axis_lab <- c(ca_amp = "Ca^{2*'+'}~amplitude~'(a.u.)'", ca_CTD50_ms = "CTD[50]~'(ms)'", ca_CTD90_ms = "CTD[90]~'(ms)'",
              ca_rise_ms = "Ca^{2*'+'}~rise~time~'(ms)'", ca_decay_ms = "Ca^{2*'+'}~decay~time~'(ms)'",
              contr_amp = "Contraction~amplitude~'(a.u.)'", contr_rise_ms = "Contraction~rise~time~'(ms)'",
              contr_decay_ms = "Contraction~decay~time~'(ms)'")
stats5 <- pmap_dfr(pairs, function(x, y) {
  d <- wm %>% select(all_of(c(x, y))) %>% drop_na()
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "spearman", exact = FALSE))
  tibble(x, y, n = nrow(d), rho = unname(ct$estimate), p = ct$p.value)
})
print(as.data.frame(stats5), digits = 3)
panels <- pmap(stats5, function(x, y, n, rho, p) {
  ggplot(wm, aes(.data[[x]], .data[[y]], fill = Condition, shape = Cell_Line)) +
    geom_point(size = 1.5, colour = "grey25", stroke = 0.25) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.4, size = 5.5 / .pt, colour = "grey20",
             parse = TRUE, label = sprintf("rho == '%.2f'*','~%s", rho, fmt_p_pm(p))) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_shape_manual(values = c("1.6" = 21, "1.6.24" = 24), labels = line_names) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.18))) +
    labs(x = parse(text = axis_lab[x]), y = parse(text = axis_lab[y])) +
    pub + theme(legend.position = "none", axis.title = element_text(size = 7))
})
leg <- cowplot::get_legend(
  ggplot(wm, aes(ca_amp, contr_amp, fill = Condition, shape = Cell_Line)) +
    geom_point(size = 2, colour = "grey25", stroke = 0.25) +
    scale_fill_manual(values = group_colors, name = "Condition",
                      guide = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_manual(values = c("1.6" = 21, "1.6.24" = 24), labels = line_names, name = "Line") +
    pub + theme(legend.position = "right", legend.text = element_text(size = 6)))
p5 <- cowplot::plot_grid(plotlist = c(panels, list(leg)), ncol = 4, align = "hv",
                         labels = c(letters[1:7], ""), label_size = 8)
ggsave(file.path(out_dir, paste0("supp_figure5_", analysis_date, ".pdf")), p5,
       width = 180, height = 95, units = "mm", device = cairo_pdf)

for (f in c("supp_figure2_", "supp_figure4_", "supp_figure5_"))
  ggsave(file.path(out_dir, paste0(f, analysis_date, ".png")),
         get(c(supp_figure2_ = "p2", supp_figure4_ = "p4", supp_figure5_ = "p5")[f]),
         width = 180, height = 95, units = "mm", dpi = 300, bg = "white")
cat("  Saved:", out_dir, "\n")
cat(sprintf("  Fig. 4 n: calcium %d, contraction %d cells. Fig. 5: %d wells.\n", c_ca$n, c_co$n, nrow(wm)))
