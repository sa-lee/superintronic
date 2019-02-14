
# this is equivalent to join_overlap_intersect()
rangeCoverage <- function(gr, cvg) {
  gr <- IRanges::reduce(gr)
  ans <- unlist(cvg[gr], use.names=FALSE)
  if (S4Vectors::runValue(BiocGenerics::strand(gr))[[1L]] == "-")
    ans <- rev(ans)
  ans
}


# can also do a PCA here

counts_by_type <- function(table, var) {
  S4Vectors::do.call(
    "cbind",
    S4Vectors::split(table[[var]], table[["source"]])
  )
}



# exon_counts <- counts_by_type(exin_by_sample, "exon.average_coverage")
# 
# pcs <- prcomp(exon_counts[1:500, ])
# 
# exon_pc_plot <- with(pcs, {
#   df <- cbind(pcs$rotation, unique(exin_by_sample[, c("Kit", "source")]))
#   ggplot2::ggplot(as.data.frame(df), 
#                   aes(x = PC1, y = PC2, colour = Kit )
#   ) +
#     ggplot2::geom_point(stroke = 1) +
#     ggplot2::coord_fixed(pcs$sdev[2] / pcs$sdev[1]) +
#     ggplot2::scale_color_brewer(palette = "Dark2") +
#     ggplot2::guides(color = FALSE)
# })
# 
# intron_counts <- counts_by_type(exin_by_sample, "intron.average_coverage")
# pcs <- prcomp(intron_counts[1:500, ])
# 
# intron_pc_plot <- with(pcs, {
#   df <- cbind(pcs$rotation, unique(exin_by_sample[, c("Kit", "source")]))
#   ggplot2::ggplot(as.data.frame(df), 
#                   aes(x = PC1, y = PC2, colour = Kit )
#   ) +
#     ggplot2::geom_point(stroke = 1) +
#     ggplot2::coord_fixed(pcs$sdev[2] / pcs$sdev[1]) +
#     ggplot2::scale_color_brewer(palette = "Dark2")
# })
# 
# pca_plots <- patchwork::wrap_plots(exon_pc_plot,
#                                    intron_pc_plot,
#                                    ncol = 2) +
#   patchwork::plot_annotation(title = "PCA of exon and intron coverage")
# 
# ggplot2::ggsave(here::here("figures", "pca.png"), pca_plots)