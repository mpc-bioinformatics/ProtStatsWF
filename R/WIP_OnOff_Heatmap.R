

#' Calculate on/off proteins
#'
#' @param D data set containg only protein intensities
#' @param id data frame containing ID columns
#' @param group factor containing the groups
#' @param max_vv_off off: < max_vv_off valid values
#' @param min_vv_on on: > min_vv_on valid values
#' @param protein_id_col column on id containing the protein IDs used for mapping
#'
#' @return data frame with number of valid values per group (absolute and relative) and on/off status
#' @export
#'
#' @examples # TODO
calculate_onoff <- function(D, id, group, max_vv_off, min_vv_on, protein_id_col = 1) {

  ### TODO: check that protein_id_col has only unique entries, otherwise the on/off calculation will fail

  group <- droplevels(group)
  nr_groups <- length(levels(group))

  #Gene.names <- id[, gene_names_col]
  Protein.IDs <- id[, protein_id_col]

  ### if gene name is empty, replace it with protein accession
  # ind_NA_gene <- which(is.na(Gene.names) | Gene.names == "")
  # if (length(ind_NA_gene) > 0) {
  #   #Gene.names[ind_NA_gene] <-
  #   id[ind_NA_gene, gene_names_col] <- Protein.IDs[ind_NA_gene]
  # }
  # id[, gene_names_col] <- make.unique(id[, gene_names_col])
  # # Gene.names <- make.unique(Gene.names)

  ## converting to long format
  D_long <- tidyr::pivot_longer(data = cbind(Protein.IDs = Protein.IDs, D), cols = colnames(D))
  D_long$group <- group[match(D_long$name, colnames(D))]

  ## calculate on/off values
  D_onoff <- D_long %>% dplyr::group_by(group, Protein.IDs) %>%
    dplyr::summarise(valid_values = sum(!is.na(value)), valid_values_rel = mean(!is.na(value)))

  ## convert to wide format again
  D_onoff_wide <- tidyr::pivot_wider(D_onoff,
                              id_cols = Protein.IDs,
                              values_from = c(valid_values, valid_values_rel),
                              names_from = group)
  ind <- match(D_onoff_wide$Protein.IDs, Protein.IDs)

  # if (relative) {
  #   cols <- colnames(D_onoff_wide)[(nr_groups+2):(2*nr_groups+1)]
  # } else {
     cols <- colnames(D_onoff_wide)[2:(nr_groups+1)]
  # }

  ### calculate, if protein is on/off
  res_onoff <- apply(D_onoff_wide[,cols], 1, function(x) {
    #if (is.null(onoffdiff))
      isonoff <- any(x <= max_vv_off) & any(x >= min_vv_on)
    #if (!is.null(onoffdiff)) isonoff <- diff(range(x)) >= onoffdiff
    return(isonoff)
  })


  D_onoff_wide_tmp <- cbind(D_onoff_wide, isonoff = res_onoff)
  D_onoff_wide_tmp2 <- cbind(id[ind,], D_onoff_wide_tmp[,-1])
  return(D_onoff_wide_tmp2)
  #write.xlsx(D_onoff_wide_tmp2, paste0(output_path, "table_onoff", suffix, ".xlsx"), keepNA = TRUE, overwrite = TRUE)

}






################################################################################
################################################################################
################################################################################
#
#
# Onoff_plus_heatmap <- function(D, id, group,
#                                gene_names_col, protein_id_col,
#                                max_vv_off, min_vv_on, onoffdiff = NULL, relative = FALSE,
#                                output_path,
#                                plot_height = 10, plot_width = 10, suffix = "", ... ){
#   require(tidyverse)
#
#
#
#   D_onoff_wide2 <- D_onoff_wide_tmp2[D_onoff_wide_tmp2$isonoff, ]
#
#   #### TODO: D_onoff_wide2 ist leer, weil isonoff für alles Falsch ist
#   ##validValue_col <- which(colnames(D_onoff_wide2) %in% cname)
#  # cname <<- colnames(D_onoff_wide2)[-c(1, ncol(D_onoff_wide2))]
#   ####
#   D_onoff2_long <- as.data.frame(pivot_longer(D_onoff_wide2, cols = cols)) # [1]:cname[length(cname)]
#   if (relative) {
#     D_onoff2_long$group <- str_replace(D_onoff2_long$name, "valid_values_rel_", "")
#   } else {
#     D_onoff2_long$group <- str_replace(D_onoff2_long$name, "valid_values_", "")
#   }
#
#   D_onoff2_long$group <- factor(D_onoff2_long$group, levels = levels(group))
#
#
#   ### TODO: clustering für Reihenfolge/order der Proteine?
#   ord <- do.call(order, args = c(as.list(D_onoff_wide2[, cols]), decreasing = TRUE))
#
#   #cl <- hclust(dist(D_onoff_wide2[, cols], method = "manhattan"), method="complete")
#   D_onoff2_long[, gene_names_col] <- factor(D_onoff2_long[, gene_names_col],
#                                      levels = D_onoff_wide2[, gene_names_col][ord])
#
#
#   D_onoff2_long_2 <<- D_onoff2_long
#
#   pl <- ggplot(data = D_onoff2_long, aes(x = group, y = Gene.names, fill = value)) +
#     geom_tile() +  ylab("Gene name") + xlab("group") + theme_bw()
#
#   #if (onoffGreaterThanEqual < 1 | !is.null(onoffdiff)) {
#   pl <- pl + scale_fill_gradient(limits = c(0,max(D_onoff2_long$value)), low = "white", high = "forestgreen") #
#   pl <- pl + theme(axis.text = element_text(size = rel(1.8)),
#                    axis.title = element_text(size = rel(1.8)),
#                    legend.title = element_text(size=rel(1.8)),
#                    legend.text = element_text(size=rel(1.8)))
#   pl
#
#
#   #print(length(unique(D_onoff2_long)))
#   #} else {
#     # colours_green <- scales::seq_gradient_pal("white", "forestgreen", "Lab")(seq(0,1,length.out=6))
#   #  colours_green <- scales::seq_gradient_pal("white", "forestgreen", "Lab")(seq(0,1,length.out=(length(unique(D_onoff2_long))+1) ))
#    # pl <- pl + scale_fill_manual(values=colours_green, name = "Valid \nvalues")
#   #}
#
#   ggsave(filename = paste0(output_path, "Heatmap_onoff", suffix,".png"),
#          device = "png", plot = pl,
#          height = plot_height, width = plot_width, dpi = 300)
# }
#
#
#
#
#
