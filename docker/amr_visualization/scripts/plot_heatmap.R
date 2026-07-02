#!/usr/bin/env Rscript
# plot_heatmap.R - R twin of plot_heatmap.py (ggtree cladogram + geom_tile heatmap
# aligned with aplot). Avoids gheatmap (broken with ggplot2 >= 3.5).
# CLI: --matrix --genes --tree --metadata --annotate --title --out --png --svg
suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(ape); library(aplot)
})

CELL <- c(resistance = "#2c5f8a", virulence = "#a93226", stress = "#ca6f1e",
          plasmid = "#148f77", other = "#7f8c8d", absent = "#ffffff")
CAT_ORDER <- c("virulence", "resistance", "stress", "plasmid", "other")
BAND_LABEL <- c(resistance = "AMR", virulence = "Virulence", stress = "Stress",
                plasmid = "Plasmid", other = "Other")
INK <- "#27333f"

ARGS <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, multiple = FALSE, default = NULL) {
  i <- which(ARGS == name); if (!length(i)) return(if (multiple) character() else default)
  v <- character(); j <- i[1] + 1
  while (j <= length(ARGS) && !startsWith(ARGS[j], "--")) { v <- c(v, ARGS[j]); j <- j + 1 }
  if (multiple) v else if (length(v)) v[1] else default
}
matrix_f <- getflag("--matrix"); genes_f <- getflag("--genes")
tree_f <- getflag("--tree"); title <- getflag("--title", default = "Resistome / virulome")
out <- getflag("--out", default = "amr_heatmap.pdf")
png_f <- getflag("--png"); svg_f <- getflag("--svg")

mat <- read.delim(matrix_f, row.names = 1, check.names = FALSE)
cat_of <- setNames(rep("resistance", ncol(mat)), colnames(mat))
if (!is.null(genes_f) && file.exists(genes_f)) {
  g <- read.delim(genes_f, check.names = FALSE)
  cat_of[as.character(g$gene)] <- as.character(g$category)
}
# order columns by category block then name
ord <- order(match(cat_of[colnames(mat)], CAT_ORDER), colnames(mat))
mat <- mat[, ord, drop = FALSE]; cats <- cat_of[colnames(mat)]

# strip labels: full name, but abbreviated to 3 chars for single-column blocks
present_cats <- intersect(CAT_ORDER, unique(cats))
ccount <- table(cats)
strip_label <- setNames(vapply(present_cats, function(c) {
  lab <- BAND_LABEL[[c]]
  if (as.integer(ccount[c]) == 1) substr(lab, 1, 3) else lab
}, character(1)), present_cats)

# tree (cladogram) -> tip order; else alphabetical
have_tree <- !is.null(tree_f) && file.exists(tree_f)
if (have_tree) {
  tr <- read.tree(tree_f)
  tp <- ggtree(tr, branch.length = "none")          # cladogram: aligned tips
  td <- tp$data[tp$data$isTip, ]
  tip_order <- td$label[order(td$y)]                 # bottom -> top
  rows <- intersect(tip_order, rownames(mat))
} else {
  rows <- sort(rownames(mat), decreasing = TRUE)
}
mat <- mat[rows, , drop = FALSE]

# long data frame: fill = category if present else "absent"
long <- data.frame(
  sample = factor(rep(rownames(mat), times = ncol(mat)), levels = rows),
  gene   = factor(rep(colnames(mat), each = nrow(mat)), levels = colnames(mat)),
  val    = as.vector(as.matrix(mat)),
  category = factor(rep(strip_label[cats], each = nrow(mat)),
                    levels = unname(strip_label[present_cats])),
  stringsAsFactors = FALSE
)
long$fill <- ifelse(long$val > 0, rep(cats, each = nrow(mat)), "absent")
long$fill <- factor(long$fill, levels = names(CELL))

ht <- ggplot(long, aes(gene, sample, fill = fill)) +
  geom_tile(color = "#eef1f3", linewidth = 0.3) +
  facet_grid(cols = vars(category), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = CELL, guide = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, colour = INK),
        axis.text.y = element_text(size = 7, colour = INK),
        panel.grid = element_blank(),
        panel.spacing = unit(3, "pt"),
        strip.text = element_text(face = "bold", colour = INK),
        strip.background = element_rect(fill = "#eef1f3", colour = NA))

w <- max(6, 2.4 + 0.22 * ncol(mat)); h <- max(3, 0.32 * nrow(mat) + 1.2)
plt <- if (have_tree) ht %>% insert_left(tp, width = 0.28) else ht

ggsave(out, plt, width = w, height = h, limitsize = FALSE)
if (!is.null(png_f)) ggsave(png_f, plt, width = w, height = h, dpi = 300, limitsize = FALSE)
if (!is.null(svg_f)) { grDevices::svg(svg_f, width = w, height = h); print(plt); dev.off() }
message(sprintf("[ok] figure -> %s", out))