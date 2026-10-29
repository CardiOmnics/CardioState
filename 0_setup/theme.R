library(ggplot2)

theme_project <- function() {
  theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 14),
      strip.background = element_blank()
    )
}

# Publication theme for journal figures, at Nature sizes: all text 5-7 pt
# (axis titles, strips and legend titles 7 pt, tick labels and legend text
# 6 pt). Panel letters are 8 pt bold lowercase (a, b, c), added with
# cowplot::plot_grid(labels = "auto", label_size = 8) or patchwork
# plot_annotation(tag_levels = "a").
journal_theme <- theme_bw(base_size = 7) +
  theme(
    axis.text        = element_text(size = 6),
    axis.title       = element_text(size = 7),
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 6),
    strip.text       = element_text(size = 7),
    strip.background = element_blank(),
    plot.title       = element_text(size = 7, face = "bold"),
    plot.subtitle    = element_text(size = 6),
    plot.caption     = element_text(size = 5.5),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.4),
    panel.border     = element_rect(linewidth = 0.4, fill = NA)
  )
