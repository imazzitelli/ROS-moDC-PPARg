# =============================================================================
# RNA-seq analysis and figures (Figures 3C-E and 4A-D, 4F-J)
# GM-CSF, GM-CSF+IL-4 and GM-CSF+NAC monocyte-derived cells (n = 4 donors)
#
# Runs end to end: kallisto -> tximport -> DESeq2 -> all RNA-seq panels.
#
# Inputs (paths relative to the project root):
#   data/kallisto/*.tsv              one kallisto abundance file per sample,
#                                    named <donor>_<condition>.tsv
#                                    (e.g. D1_GM.tsv, D1_GM+IL4.tsv, D1_GM+NAC.tsv)
#   data/transcripts_to_genes.csv    transcript-to-gene table for tximport
#
# Outputs:
#   figures/   one PDF per panel
#   results/   DESeq2 tables, TF activity and GSEA results, sessionInfo
#
# Note: decoupleR::get_collectri() and msigdbr download/ship versioned
# databases; results can shift slightly with other versions (see sessionInfo).
# =============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(viridis)
  library(factoextra)
  library(matrixStats)
  library(ComplexHeatmap)
  library(circlize)
  library(decoupleR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Colors and labels used across the paper
pal.cond <- c("GM" = "#B3B3B3", "GM_IL4" = "#A9DCF2", "GM_NAC" = "#5FBFAB")
lab.cond <- c("GM" = "GM-CSF", "GM_IL4" = "GM-CSF+IL-4", "GM_NAC" = "GM-CSF+NAC")
pal.donor <- setNames(viridis(4, option = "G", begin = .25, end = .90),
                      c("D1", "D2", "D3", "D4"))


# =============================================================================
# 1. Import (kallisto + tximport)
# =============================================================================

files <- list.files("data/kallisto", pattern = "\\.tsv$", full.names = TRUE)
names(files) <- basename(files)

# Sample metadata is parsed from the file names, not from their position:
# list.files() sorts according to the system locale, so the file order
# differs between Windows, macOS and Linux.
cond.raw <- sub("^D[0-9]+_(.*)\\.tsv$", "\\1", names(files))
Metadata <- data.frame(
  samples.name = names(files),
  donor        = sub("^(D[0-9]+)_.*$", "\\1", names(files)),
  experimental.condition = dplyr::recode(cond.raw,
                                         "GM"     = "GM",
                                         "GM+IL4" = "GM_IL4",
                                         "GM+NAC" = "GM_NAC",
                                         .default = "Monocytes")
)

# Canonical order: donor, then GM, GM_IL4, GM_NAC, Monocytes
ord <- order(Metadata$donor,
             factor(Metadata$experimental.condition,
                    levels = c("GM", "GM_IL4", "GM_NAC", "Monocytes")))
Metadata <- Metadata[ord, ]
files    <- files[Metadata$samples.name]
rownames(Metadata) <- Metadata$samples.name

stopifnot(nrow(Metadata) == 16,
          all(table(Metadata$experimental.condition) == 4),
          all(table(Metadata$donor) == 4))
write.csv(Metadata, "results/Metadata.csv", row.names = FALSE)

tx2gene <- read.csv("data/transcripts_to_genes.csv")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                ignoreAfterBar = TRUE)
write.csv(txi$counts, "results/Count.Matrix.csv")


# =============================================================================
# 2. DESeq2
# =============================================================================

coldata <- Metadata
coldata$experimental.condition <- factor(coldata$experimental.condition)
coldata$donor <- factor(coldata$donor)

dds <- DESeqDataSetFromMatrix(countData = round(txi$counts),
                              colData   = coldata,
                              design    = ~ donor + experimental.condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # genes with >= 10 reads in total
dds <- DESeq(dds)
saveRDS(dds, "results/dds.rds")

norm.counts <- counts(dds, normalized = TRUE)
write.csv(norm.counts, "results/Normalized.count.matrix.DESeq2.csv")

mat.vst <- assay(vst(dds, blind = FALSE))

# Library sizes (QC, not in the paper)
qc <- data.frame(sample = colnames(dds),
                 raw = colSums(counts(dds)),
                 normalized = colSums(norm.counts))
write.csv(qc, "results/QC_library_sizes.csv", row.names = FALSE)

# Differential expression (FDR < 0.05, |log2FC| > 0.58)
padj.cutoff <- 0.05
lfc.cutoff  <- 0.58

run_contrast <- function(test, ref) {
  res <- results(dds, contrast = c("experimental.condition", test, ref))
  tab <- as.data.frame(res[order(res$pvalue), ])
  deg <- subset(tab, !is.na(padj) & padj < padj.cutoff &
                  abs(log2FoldChange) > lfc.cutoff)
  write.csv(tab, sprintf("results/DESeq2_%s_vs_%s.csv", test, ref))
  write.csv(deg, sprintf("results/DEG_%s_vs_%s.csv",    test, ref))
  list(res = res, tab = tab, deg = deg)
}

NAC.vs.IL4 <- run_contrast("GM_NAC", "GM_IL4")
NAC.vs.GM  <- run_contrast("GM_NAC", "GM")
IL4.vs.GM  <- run_contrast("GM_IL4", "GM")


# =============================================================================
# 3. Donor-wise z-scores (used by Figure 4A and 4B)
#    For each donor, each gene is z-scored across that donor's 4 samples
#    (including Monocytes). Genes with zero variance in a donor give NaN.
# =============================================================================

z.all <- do.call(cbind, lapply(c("D1", "D2", "D3", "D4"), function(d) {
  idx <- Metadata$donor == d
  t(scale(t(norm.counts[, idx, drop = FALSE])))
}))
stopifnot(identical(colnames(z.all), Metadata$samples.name))


# =============================================================================
# FIGURE 3C-D. Cytokines and chemokines, GM-CSF+NAC vs GM-CSF+IL-4
# FIGURE 3E.   Cytokines and chemokines, GM-CSF+NAC vs GM-CSF
# =============================================================================

# Chemokines: defined by nomenclature
chemo <- grep("^(CCL|CXCL|XCL|CX3CL)[0-9]+$|^(PPBP|PF4|PF4V1)$",
              rownames(mat.vst), value = TRUE)

# Cytokines: GO:0005125 "cytokine activity" (including child terms),
# without chemokines
cyto <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0005125",
                              keytype = "GOALL", columns = "SYMBOL")$SYMBOL
cyto <- setdiff(intersect(unique(cyto), rownames(mat.vst)), chemo)

# Expressed genes: mean normalized counts >= 20
expr.ok <- rownames(norm.counts)[rowMeans(norm.counts) >= 20]
chemo   <- intersect(chemo, expr.ok)
cyto    <- intersect(cyto,  expr.ok)

# Annotated under GO:0005125 but primarily enzymes, protease inhibitors
# or morphogens; excluded from the cytokine panels
cyto <- setdiff(cyto, c("NAMPT", "EDN1", "TIMP1", "CMTM1", "WNT5A", "WNT5B"))

# Lollipop: genes with padj < 0.05 and |log2FC| > 1
make_lollipop <- function(res, genes, title, base_size = 15, font_y = 13) {
  df <- as.data.frame(res[intersect(genes, rownames(res)), ])
  df$gene <- rownames(df)
  df <- subset(df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
  df <- df[order(df$log2FoldChange), ]
  df$gene <- factor(df$gene, levels = df$gene)

  ggplot(df, aes(x = log2FoldChange, y = gene)) +
    geom_segment(aes(x = 0, xend = log2FoldChange, yend = gene),
                 color = "grey70", linewidth = .4) +
    geom_point(aes(fill = log2FoldChange, size = -log10(padj)),
               shape = 21, color = "black", stroke = .3) +
    scale_fill_viridis_c(name = "log2FC") +
    scale_size_continuous(name = "-log10(padj)", range = c(2, 6)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    labs(title = title, x = "log2FC", y = NULL) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.major.y = element_line(linewidth = .2),
          panel.grid.minor   = element_blank(),
          axis.text.y        = element_text(face = "italic", size = font_y))
}

ggsave("figures/Fig3C_cytokines_NACvsIL4.pdf",
       make_lollipop(NAC.vs.IL4$res, cyto,
                     "Cytokines - GM-CSF+NAC vs GM-CSF+IL4"),
       width = 7, height = 8)

ggsave("figures/Fig3D_chemokines_NACvsIL4.pdf",
       make_lollipop(NAC.vs.IL4$res, chemo,
                     "Chemokines - GM-CSF+NAC vs GM-CSF+IL4"),
       width = 7, height = 6.5)

ggsave("figures/Fig3E_cytokines_NACvsGM.pdf",
       make_lollipop(NAC.vs.GM$res, cyto,
                     "Cytokines - GM-CSF+NAC vs GM-CSF"),
       width = 5.5, height = 8)

ggsave("figures/Fig3E_chemokines_NACvsGM.pdf",
       make_lollipop(NAC.vs.GM$res, chemo,
                     "Chemokines - GM-CSF+NAC vs GM-CSF",
                     base_size = 11, font_y = 8),
       width = 7, height = 6.5)


# =============================================================================
# FIGURE 4A. PCA on donor-wise z-scores (GM, GM+IL4, GM+NAC)
# =============================================================================

keep.act <- Metadata$experimental.condition != "Monocytes"
z.pca    <- na.omit(z.all)[, keep.act]          # genes finite in all donors
pca      <- prcomp(t(z.pca), center = TRUE, scale. = FALSE)
var.exp  <- 100 * pca$sdev^2 / sum(pca$sdev^2)

cond.pca <- factor(Metadata$experimental.condition[keep.act],
                   levels = c("GM", "GM_IL4", "GM_NAC"))

p.pca <- fviz_pca_ind(pca, geom.ind = "point", pointshape = 21,
                      pointsize = 4, stroke = 0.6,
                      fill.ind = cond.pca, col.ind = "black",
                      palette = unname(pal.cond), addEllipses = FALSE,
                      label = "none", axes.linetype = NA) +
  scale_fill_manual(values = unname(pal.cond), labels = unname(lab.cond),
                    name = NULL) +
  labs(title = NULL,
       x = sprintf("PC1 (%.1f%%)", var.exp[1]),
       y = sprintf("PC2 (%.1f%%)", var.exp[2])) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

ggsave("figures/Fig4A_PCA.pdf", p.pca, width = 5.5, height = 5)


# =============================================================================
# FIGURE 4B. Top 350 variable genes, GM-CSF+IL-4 vs GM-CSF+NAC
#   z-scores computed within donor across all 4 conditions, subset to the
#   8 IL-4/NAC samples and re-centered per gene. Rows and columns are
#   clustered as in gplots::heatmap.2 (Euclidean distance, complete linkage,
#   dendrograms reordered by row/column means).
# =============================================================================

sel.2 <- Metadata$experimental.condition %in% c("GM_IL4", "GM_NAC")
z.2   <- z.all[, sel.2]
top350 <- head(order(rowVars(z.2), decreasing = TRUE), 350)
m.b    <- z.2[top350, ]
m.b    <- m.b - rowMeans(m.b)
meta.b <- droplevels(Metadata[sel.2, ])

row.dend <- reorder(as.dendrogram(hclust(dist(m.b))), rowMeans(m.b))
col.dend <- reorder(as.dendrogram(hclust(dist(t(m.b)))), colMeans(m.b))

zlim  <- max(abs(m.b))
hm.b <- Heatmap(
  m.b, name = "z-score",
  col = colorRamp2(seq(-zlim, zlim, length.out = 50), viridis(50)),
  cluster_rows = rev(row.dend), cluster_columns = col.dend,
  show_column_dend = FALSE, row_dend_width = unit(12, "mm"),
  show_row_names = FALSE, show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(
    Condition = factor(meta.b$experimental.condition),
    col = list(Condition = pal.cond[c("GM_IL4", "GM_NAC")])),
  bottom_annotation = HeatmapAnnotation(
    Donor = meta.b$donor, col = list(Donor = pal.donor)),
  column_title = "Top 350 variable genes\nGM-CSF+IL4 vs GM-CSF+NAC",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE)

pdf("figures/Fig4B_heatmap_top350_IL4_vs_NAC.pdf", width = 6, height = 8)
draw(hm.b, merge_legend = TRUE)
dev.off()


# =============================================================================
# FIGURE 4C. Volcano plot, GM-CSF+NAC vs GM-CSF+IL-4
# =============================================================================

plot_volcano <- function(df, title, pCutoff = 0.05, FCcutoff = 0.58,
                         xlim_val = c(-10, 10), top_n_labels = 15) {
  pal <- viridis(3)
  color_map <- c("Downregulated" = pal[2], "NS" = "grey80", "Upregulated" = pal[1])

  df <- df %>%
    mutate(gene = rownames(df),
           negLog10P = -log10(padj),
           category = case_when(
             padj < pCutoff & log2FoldChange >  FCcutoff ~ "Upregulated",
             padj < pCutoff & log2FoldChange < -FCcutoff ~ "Downregulated",
             TRUE ~ "NS"),
           category = factor(category,
                             levels = c("Downregulated", "NS", "Upregulated")))

  labels_df <- df %>%
    filter(category != "NS") %>%
    arrange(padj) %>%
    group_by(category) %>%
    slice_head(n = top_n_labels) %>%
    ungroup()

  n_up   <- sum(df$category == "Upregulated")
  n_down <- sum(df$category == "Downregulated")
  cat(title, "-> up:", n_up, "| down:", n_down, "\n")

  ggplot(df, aes(x = log2FoldChange, y = negLog10P, color = category)) +
    geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(pCutoff), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_point(size = 1.6, alpha = 0.85) +
    geom_label_repel(data = labels_df, aes(label = gene), size = 2.6,
                     max.overlaps = 100, segment.size = 0.2,
                     segment.color = "grey30", label.size = 0.25,
                     label.padding = 0.15, fill = "white", box.padding = 0.3,
                     show.legend = FALSE, seed = 42) +
    scale_color_manual(values = color_map, name = NULL) +
    coord_cartesian(xlim = xlim_val) +
    labs(title = title,
         x = expression(Log[2]~"Fold Change"),
         y = expression(-Log[10]~"FDR"),
         caption = sprintf("%d up / %d down  (FDR < %s, |log2FC| > %s)",
                           n_up, n_down, pCutoff, FCcutoff)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top",
          legend.text = element_text(size = 9),
          plot.caption = element_text(hjust = 0.5, size = 9, face = "italic"),
          panel.grid.minor = element_blank())
}

ggsave("figures/Fig4C_volcano_NACvsIL4.pdf",
       plot_volcano(NAC.vs.IL4$tab, "GM-CSF+NAC vs GM-CSF+IL-4"),
       width = 6, height = 5.5)


# =============================================================================
# FIGURE 4D. Transcription factor activity (decoupleR ULM, CollecTRI)
#   ULM on vst data without donor correction; paired test on donor-matched
#   NAC - IL4 differences; donor effect removed only for display.
# =============================================================================

meta.act <- Metadata[keep.act, ]
meta.act$experimental.condition <- factor(meta.act$experimental.condition,
                                          levels = c("GM", "GM_IL4", "GM_NAC"))
meta.act$donor <- factor(meta.act$donor)

net  <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
acts <- decoupleR::run_ulm(mat = mat.vst[, keep.act, drop = FALSE], net = net,
                           .source = "source", .target = "target",
                           .mor = "mor", minsize = 5)

act.mat <- acts %>%
  dplyr::filter(statistic == "ulm") %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
act.mat <- act.mat[, rownames(meta.act), drop = FALSE]
stopifnot(identical(colnames(act.mat), rownames(meta.act)))

# Donor-matched NAC - IL4 differences
sel.tf <- meta.act$experimental.condition %in% c("GM_IL4", "GM_NAC")
act2   <- act.mat[, sel.tf, drop = FALSE]
meta2  <- droplevels(meta.act[sel.tf, ])
donors <- levels(meta2$donor)

id.nac <- vapply(donors, function(d) which(meta2$donor == d &
                   meta2$experimental.condition == "GM_NAC"), integer(1))
id.il4 <- vapply(donors, function(d) which(meta2$donor == d &
                   meta2$experimental.condition == "GM_IL4"), integer(1))

D <- act2[, id.nac, drop = FALSE] - act2[, id.il4, drop = FALSE]  # + = higher in NAC
res.tf <- data.frame(TF   = rownames(D),
                     diff = rowMeans(D),
                     p    = apply(D, 1, function(x) t.test(x)$p.value),
                     row.names = NULL)
res.tf$padj <- p.adjust(res.tf$p, method = "BH")
res.tf <- res.tf[order(res.tf$p), ]
write.csv(res.tf, "results/TF_activity_NAC_vs_IL4_paired.csv", row.names = FALSE)

# Top 15 significant TFs per direction
sig    <- subset(res.tf, !is.na(padj) & padj < 0.05)
up.nac <- head(sig$TF[sig$diff > 0][order(-sig$diff[sig$diff > 0])], 15)
up.il4 <- head(sig$TF[sig$diff < 0][order( sig$diff[sig$diff < 0])], 15)

# Display filter: drop TFs not expressed (mean normalized counts <= 10);
# AP1 and NFKB are complexes, not genes
expressed <- rownames(norm.counts)[rowMeans(norm.counts[, keep.act]) > 10]
genes.tf  <- c(up.nac, up.il4)
keep.tf   <- genes.tf %in% c(expressed, "AP1", "NFKB")
split.tf  <- factor(rep(c("Up in NAC", "Up in IL4"),
                        c(length(up.nac), length(up.il4))),
                    levels = c("Up in NAC", "Up in IL4"))[keep.tf]
genes.tf  <- genes.tf[keep.tf]

act.disp <- limma::removeBatchEffect(act2, batch = meta2$donor)
mt  <- act.disp[genes.tf, , drop = FALSE]
mt  <- mt - rowMeans(mt)
lim <- as.numeric(quantile(abs(mt), 0.98))

stars <- cut(res.tf$padj[match(genes.tf, res.tf$TF)],
             c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", ""))

hm.tf <- Heatmap(
  mt, name = "TF activity\n(centered)",
  col = colorRamp2(seq(-lim, lim, length.out = 11), viridis(11)),
  top_annotation = HeatmapAnnotation(
    Condition = meta2$experimental.condition,
    Donor     = meta2$donor,
    col = list(Condition = pal.cond[c("GM_IL4", "GM_NAC")], Donor = pal.donor),
    annotation_name_gp = gpar(fontsize = 10)),
  right_annotation = rowAnnotation(
    FDR = anno_text(as.character(stars), gp = gpar(fontsize = 10)),
    annotation_name_gp = gpar(fontsize = 9)),
  row_split = split.tf, column_split = meta2$experimental.condition,
  border = TRUE, border_gp = gpar(col = "black", lwd = 1.2),
  row_gap = unit(1.5, "mm"), column_gap = unit(1.5, "mm"),
  cluster_rows = FALSE, cluster_row_slices = FALSE,
  cluster_columns = TRUE, cluster_column_slices = FALSE,
  clustering_distance_columns = "pearson",
  clustering_method_columns   = "ward.D2",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  row_title_gp = gpar(fontsize = 11),
  column_title_gp = gpar(fontsize = 11),
  show_column_names = FALSE)

pdf("figures/Fig4D_TF_activity_NACvsIL4.pdf", width = 6.8, height = 8)
draw(hm.tf, merge_legend = TRUE,
     column_title = "TF activity - GM-CSF+NAC vs GM-CSF+IL-4",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"))
dev.off()


# =============================================================================
# FIGURE 4F and 4H. Preranked GSEA on MSigDB Hallmark
#   Genes with mean normalized counts >= 50 in the two conditions compared;
#   log2FC and paired t-test on log2(normalized counts + 1), ranked by
#   sign(log2FC) * -log10(p).
# =============================================================================

hallmark.df <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                        error = function(e) msigdbr(species = "Homo sapiens",
                                                    category = "H"))
hallmark.sets <- split(hallmark.df$gene_symbol, hallmark.df$gs_name)

cols_of <- function(cond) {
  s <- Metadata[Metadata$experimental.condition == cond, ]
  s$samples.name[order(s$donor)]            # donor-matched order
}

pretty_term <- function(x) {
  s <- str_to_title(gsub("_", " ", sub("^HALLMARK_", "", x)))
  for (ac in c("DNA", "RNA", "ROS", "TGF", "IL2", "E2F", "G2M", "KRAS", "UV"))
    s <- gsub(paste0("\\b", str_to_title(ac), "\\b"), ac, s)
  s
}

run_gsea <- function(test, ref, file_tag) {
  M <- norm.counts[, c(cols_of(test), cols_of(ref))]
  meanexpr <- rowMeans(M)
  M <- M[meanexpr >= 50, ]

  log_t <- log2(M[, cols_of(test)] + 1)
  log_r <- log2(M[, cols_of(ref)]  + 1)

  de <- data.frame(
    gene   = rownames(M),
    log2FC = rowMeans(log_t - log_r),
    p      = vapply(seq_len(nrow(M)), function(i) {
      x <- log_t[i, ]; y <- log_r[i, ]
      if (sd(x - y) == 0) return(NA_real_)
      tryCatch(t.test(x, y, paired = TRUE)$p.value, error = function(e) NA_real_)
    }, numeric(1)))
  de$padj <- p.adjust(de$p, method = "BH")
  write.csv(de, sprintf("results/paired_DE_%s.csv", file_tag), row.names = FALSE)

  rnk  <- subset(de, is.finite(log2FC) & is.finite(p) & p > 0)
  stat <- sort(setNames(sign(rnk$log2FC) * -log10(rnk$p), rnk$gene),
               decreasing = TRUE)

  set.seed(42)
  fg <- fgseaMultilevel(pathways = hallmark.sets, stats = stat,
                        minSize = 10, maxSize = 500, eps = 0)
  fg <- fg[order(fg$NES), ]

  out <- as.data.frame(fg)
  out$Term <- pretty_term(out$pathway)
  out$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", character(1))
  write.csv(out[, c("pathway", "Term", "ES", "NES", "pval", "padj", "size", "leadingEdge")],
            sprintf("results/GSEA_Hallmark_%s.csv", file_tag), row.names = FALSE)
  out
}

plot_gsea <- function(gs, title, x_lab, ktop = 20) {
  gs  <- subset(gs, padj < 0.25)
  up  <- head(gs[order(-gs$NES), ], ktop)
  dn  <- head(gs[order( gs$NES), ], ktop)
  d   <- bind_rows(dn, up) %>%
    distinct(Term, .keep_all = TRUE) %>%
    arrange(NES) %>%
    mutate(Term     = factor(Term, levels = Term),
           label    = sprintf("%+.2f%s", NES, ifelse(padj < 0.05, "*", "")),
           fontface = ifelse(padj < 0.05, "bold", "plain"))

  p <- ggplot(d, aes(x = NES, y = Term, fill = NES)) +
    geom_col(width = 0.74, colour = "black", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.5) +
    geom_text(aes(label = label, fontface = fontface,
                  hjust = ifelse(NES < 0, 1.15, -0.15)),
              size = 2.9, colour = "grey10") +
    scale_fill_viridis_c(option = "D", name = "NES") +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
    labs(x = x_lab, y = NULL, title = title,
         subtitle = "Top pathways at FDR < 0.25; * = FDR < 0.05") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 12, colour = "grey30"),
          axis.text.y = element_text(size = 12),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4))
  list(plot = p, n = nrow(d))
}

gsea.NAC.IL4 <- run_gsea("GM_NAC", "GM_IL4", "NAC_vs_IL4")
p4F <- plot_gsea(gsea.NAC.IL4,
                 "GM-CSF+NAC vs GM-CSF+IL4 \u2014 GSEA (Hallmark Pathways)",
                 "Normalized Enrichment Score (NES)\n\u2190 enriched in GM-CSF+IL4       enriched in GM-CSF+NAC \u2192")
ggsave("figures/Fig4F_GSEA_Hallmark_NACvsIL4.pdf", p4F$plot,
       width = 8, height = 7, device = cairo_pdf)

gsea.NAC.GM <- run_gsea("GM_NAC", "GM", "NAC_vs_GM")
p4H <- plot_gsea(gsea.NAC.GM,
                 "GM-CSF vs GM-CSF+NAC \u2014 GSEA (Hallmark Pathways)",
                 "Normalized Enrichment Score (NES)\n\u2190 enriched in GM-CSF      enriched in GM-CSF+NAC \u2192")
ggsave("figures/Fig4H_GSEA_Hallmark_NACvsGM.pdf", p4H$plot,
       width = 9.6, height = 0.42 * p4H$n + 1.8, device = cairo_pdf)


# =============================================================================
# FIGURE 4G. Hallmark over-representation, genes up in GM-CSF+NAC vs IL-4
#   All DEGs (FDR < 0.05, |log2FC| > 0.58) with log2FC > 0.5; background =
#   all genes tested by DESeq2 (non-NA padj).
# =============================================================================

col.entrez <- intersect(c("ncbi_gene", "entrez_gene"), colnames(hallmark.df))[1]
m_t2g <- as.data.frame(hallmark.df[, c("gs_name", col.entrez)])
colnames(m_t2g) <- c("gs_name", "entrez_gene")

up.genes <- rownames(subset(NAC.vs.IL4$deg, log2FoldChange > 0.5))
ez.up    <- unique(bitr(up.genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
universe <- unique(bitr(rownames(subset(NAC.vs.IL4$tab, !is.na(padj))),
                        "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)

ora.hallmark <- enricher(ez.up, universe = universe, TERM2GENE = m_t2g,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
write.csv(as.data.frame(ora.hallmark), "results/ORA_Hallmark_up_NAC_vs_IL4.csv",
          row.names = FALSE)

df.ora <- as.data.frame(ora.hallmark) %>%
  filter(p.adjust < 0.05) %>%
  slice_min(p.adjust, n = 20) %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))),
         Description   = gsub("_", " ", gsub("^HALLMARK_", "", Description)),
         Description   = factor(Description,
                                levels = Description[order(GeneRatio_num)]))

p4G <- ggplot(df.ora, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(fill = p.adjust, size = Count), shape = 21,
             color = "black", stroke = 0.4) +
  scale_fill_viridis_c(name = "p.adjust", option = "D", direction = -1) +
  scale_size_continuous(name = "Count", range = c(4, 9)) +
  labs(title = "Hallmark Pathways \u2013 Upregulated in GM-CSF+NAC",
       x = "GeneRatio", y = NULL) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 12, margin = margin(t = 8)),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank())

ggsave("figures/Fig4G_ORA_Hallmark_up_NAC.pdf", p4G,
       width = 8, height = 6, device = cairo_pdf)


# =============================================================================
# FIGURE 4I. NRF2-associated genes, GM-CSF+NAC and GM-CSF+IL-4 vs GM-CSF
#   Donor-matched log2FC on DESeq2-normalized counts; bars = mean, error = SEM
# =============================================================================

nrf2 <- c("GLRX", "HMOX1", "GCLM", "NQO1", "TXN", "FTL",
          "TXNRD1", "GSR", "PRDX1", "GPX4")
nrf2 <- intersect(nrf2, rownames(norm.counts))

log_gm <- log2(norm.counts[nrf2, cols_of("GM")] + 1)
lfc_vs_gm <- function(cond) {
  m <- log2(norm.counts[nrf2, cols_of(cond)] + 1) - log_gm
  colnames(m) <- c("D1", "D2", "D3", "D4")
  as.data.frame(m) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "donor", values_to = "lfc") %>%
    mutate(cond = cond)
}
fc.nrf2 <- bind_rows(lfc_vs_gm("GM_NAC"), lfc_vs_gm("GM_IL4"))

ord.nrf2 <- fc.nrf2 %>%
  filter(cond == "GM_NAC") %>%
  group_by(gene) %>%
  summarise(m = mean(lfc), .groups = "drop") %>%
  arrange(m) %>%
  pull(gene)
fc.nrf2$gene <- factor(fc.nrf2$gene, levels = rev(ord.nrf2))
fc.nrf2$cond <- factor(fc.nrf2$cond, levels = c("GM_NAC", "GM_IL4"))

p4I <- ggplot(fc.nrf2, aes(x = lfc, y = gene, fill = cond)) +
  stat_summary(fun = mean, geom = "col", position = position_dodge(.72),
               width = .68, color = "grey20", linewidth = .3) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(.72), width = .22,
               linewidth = .4, color = "grey20") +
  geom_vline(xintercept = 0, linewidth = .5, color = "grey20") +
  scale_fill_manual(values = pal.cond[c("GM_NAC", "GM_IL4")],
                    labels = c("GM-CSF+NAC", "GM-CSF+IL-4"), name = NULL) +
  scale_x_continuous(breaks = pretty_breaks(6),
                     expand = expansion(mult = c(.05, .05))) +
  labs(x = expression(log[2]*" fold change vs GM-CSF"), y = NULL,
       title = "NRF2-associated transcriptional program") +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(face = "italic", size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, color = "grey88"),
        axis.line.x = element_line(linewidth = .4, color = "grey30"),
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        legend.position = "top",
        legend.justification = "left",
        legend.key.size = unit(4, "mm"))

ggsave("figures/Fig4I_NRF2_program.pdf", p4I, width = 6.5, height = 6.5,
       device = cairo_pdf)


# =============================================================================
# FIGURE 4J. Cholesterol/sterol genes, GM-CSF+NAC vs GM-CSF (DESeq2 log2FC)
# =============================================================================

genes_chol <- tribble(
  ~gene,     ~process,
  "MVK",     "Biosynthesis - mevalonate",
  "IDI1",    "Biosynthesis - mevalonate",
  "ACSS2",   "Biosynthesis - mevalonate",
  "HMGCR",   "Biosynthesis - mevalonate",
  "HMGCS1",  "Biosynthesis - mevalonate",
  "PMVK",    "Biosynthesis - mevalonate",
  "FDPS",    "Biosynthesis - mevalonate",
  "FDFT1",   "Biosynthesis - mevalonate",
  "MVD",     "Biosynthesis - mevalonate",
  "SC5D",    "Biosynthesis - sterols",
  "CYP51A1", "Biosynthesis - sterols",
  "HSD17B7", "Biosynthesis - sterols",
  "EBP",     "Biosynthesis - sterols",
  "SQLE",    "Biosynthesis - sterols",
  "LSS",     "Biosynthesis - sterols",
  "NSDHL",   "Biosynthesis - sterols",
  "LDLR",    "Uptake (import)",
  "STX5",    "Transport/efflux",
  "STARD4",  "Transport/efflux",
  "ABCA2",   "Transport/efflux",
  "TMEM97",  "Regulation (sensing)",
  "SREBF2",  "Regulation (sensing)"
)
process_order <- c("Biosynthesis - mevalonate", "Biosynthesis - sterols",
                   "Uptake (import)", "Transport/efflux", "Regulation (sensing)")

df.chol <- genes_chol %>%
  left_join(data.frame(gene = rownames(NAC.vs.GM$tab),
                       log2FC = NAC.vs.GM$tab$log2FoldChange), by = "gene") %>%
  filter(!is.na(log2FC))

proc_stats <- df.chol %>%
  group_by(process) %>%
  summarise(n = n(), med = median(log2FC), .groups = "drop") %>%
  mutate(leg_label = sprintf("%s  (n=%d, med %.2f)", process, n, med))

df.chol <- df.chol %>%
  left_join(proc_stats, by = "process") %>%
  mutate(process   = factor(process, levels = process_order),
         leg_label = factor(leg_label,
                            levels = proc_stats$leg_label[match(process_order,
                                                                proc_stats$process)])) %>%
  arrange(process, desc(log2FC)) %>%
  mutate(y_pos = rev(row_number()),
         lab_x = ifelse(log2FC < 0, log2FC - 0.02, log2FC + 0.02),
         lab_h = ifelse(log2FC < 0, 1, 0))

xmin <- min(df.chol$log2FC) - 0.30
xmax <- max(df.chol$log2FC) + 0.25

rects <- df.chol %>%
  group_by(process) %>%
  summarise(ymin = min(y_pos) - 0.5, ymax = max(y_pos) + 0.5, .groups = "drop") %>%
  mutate(fill_bg = rep(c("grey96", "white"), length.out = n()))
y_lines <- sort(unique(c(rects$ymin, rects$ymax)))
y_lines <- y_lines[-c(1, length(y_lines))]

p4J <- ggplot(df.chol) +
  geom_rect(data = rects, aes(ymin = ymin, ymax = ymax, xmin = xmin,
                              xmax = xmax, fill = fill_bg),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_hline(yintercept = y_lines, linetype = "dotted",
             color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = log2FC, y = y_pos, yend = y_pos,
                   color = leg_label), linewidth = 4.5) +
  geom_text(aes(x = lab_x, y = y_pos, label = gene, hjust = lab_h),
            size = 3.2, fontface = "bold", color = "grey15") +
  scale_color_viridis_d(option = "D", end = 0.85,
                        name = "Process (cholesterol/sterol)") +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.8, 0.8))) +
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0)) +
  labs(title = "Cholesterol/sterol genes - GM-CSF+NAC vs GM-CSF",
       x = expression(log[2]*" fold change (GM-CSF+NAC / GM-CSF)"), y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 8.5))

ggsave("figures/Fig4J_cholesterol_genes_NACvsGM.pdf", p4J,
       width = 7, height = 8, device = cairo_pdf)


# =============================================================================
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("Done. Figures in figures/, tables in results/\n")
