

#' Calculate on/off proteins
#'
#' @param D                \strong{data.frame} \cr
#'                         The data set containing only protein intensities of the sample.
#' @param id               \strong{data.frame} \cr
#'                         The corresponding ID columns for the parameter D.
#' @param group            \strong{character factor} \cr
#'                         The groups of the data.
#' @param max_vv_off       \strong{integer} \cr
#'                         The proteins below the threshold of the maximum number of valid values are considered off proteins.
#' @param min_vv_on        \strong{integer} \cr
#'                         The proteins above the threshold of the minimum number of valid values are considered on proteins.
#' @param protein_id_col   \strong{integer} \cr
#'                         The column in id parameter containing the protein IDs used for mapping.
#'
#' @return A data.frame with the number of valid values per group (absolute and relative) and on/off status
#'
#'
# @seealso [Onoff_plus_heatmap()]
#'
#' @examples
#'

calculate_onoff <- function(D, id, group, max_vv_off, min_vv_on, protein_id_col = 1) {

  ### TODO: check that protein_id_col has only unique entries, otherwise the on/off calculation will fail

  group <- droplevels(group)
  nr_groups <- length(levels(group))

  #Gene.names <- id[, gene_names_col]
  Protein.IDs <- id[, protein_id_col]

  ## converting to long format
  D_long <- tidyr::pivot_longer(data = cbind(Protein.IDs = Protein.IDs, D), cols = colnames(D))
  D_long$group <- group[match(D_long$name, colnames(D))]


  value <- valid_values <- valid_values_rel <- NULL # silence notes when checking the package

  ## calculate on/off values
  D_onoff <- D_long %>% dplyr::group_by(group, Protein.IDs) %>%
    dplyr::summarise(valid_values = sum(!is.na(value)), valid_values_rel = mean(!is.na(value)))

  ## convert to wide format again
  D_onoff_wide <- tidyr::pivot_wider(D_onoff,
                                     id_cols = Protein.IDs,
                                     values_from = c(valid_values, valid_values_rel),
                                     names_from = group)
  ind <- match(D_onoff_wide$Protein.IDs, Protein.IDs)

  cols <- colnames(D_onoff_wide)[2:(nr_groups+1)]


  ### calculate, if protein is on/off
  res_onoff <- apply(D_onoff_wide[,cols], 1, function(x) {
    isonoff <- any(x <= max_vv_off) & any(x >= min_vv_on)
    return(isonoff)
  })


  D_onoff_wide_tmp <- cbind(D_onoff_wide, isonoff = res_onoff)
  D_onoff_wide_tmp2 <- cbind(id[ind,], D_onoff_wide_tmp[,-1])
  return(D_onoff_wide_tmp2)
}






################################################################################
################################################################################
################################################################################
#

# Calculate on/off proteins
#
# @param RES_onoff             \strong{?} \cr
#                              ?
# @param protein_name_column   \strong{character} \cr
#                              The column name of the proteins.
# @param relative              \strong{logical} \cr
#                              If \code{TRUE}, ?
#
# @return A data.frame with the number of valid values per group (absolute and relative) and on/off status
# @export
#
# @seealso [calculate_onoff()]
#
# @examples
#

Onoff_plus_heatmap <- function(RES_onoff,
                               protein_name_column = "Gene.names",
                               relative = FALSE){

  require(tidyverse)

  ## choose only the rows with on/off proteins
  RES_onoff2 <- RES_onoff[RES_onoff$isonoff, ]

  #### TODO: D_onoff_wide2 ist leer, weil isonoff für alles Falsch ist

  validvalue_cols <- setdiff(colnames(RES_onoff2)[grep("valid_values_", colnames(RES_onoff2))], colnames(RES_onoff2)[grep("valid_values_rel_", colnames(RES_onoff2))])

  ### TODO: schlauere Methode um doppelte Proteinname zu behandeln? (Das sind meist die leeren! -> überschreiben mit protein accession z.B.)
  RES_onoff2[, protein_name_column] <- make.names(RES_onoff2[, protein_name_column], unique = TRUE)


<<<<<<< HEAD
  RES_onoff2_long <- as.data.frame(tidyr::pivot_longer(RES_onoff2, cols = tidyselect::all_of(validvalue_cols), names_to = "group"))
=======
  RES_onoff2_long <- as.data.frame(tidyr::pivot_longer(RES_onoff2, cols = tidyr::all_of(validvalue_cols), names_to = "group"))
>>>>>>> cdea481 (clean up to get rid of most warnings and notes when checking)


  if (relative) {
    RES_onoff2_long$group <- stringr::str_replace(RES_onoff2_long$group, "valid_values_rel_", "")
  } else {
    RES_onoff2_long$group <- stringr::str_replace(RES_onoff2_long$group, "valid_values_", "")
  }

  ### TODO: level Reihenfolge der Gruppe nutzen statt alphabetisch
  #RES_onoff2_long$group <- factor(RES_onoff2_long$group, levels = levels(group))


  ### TODO: clustering für Reihenfolge/order der Proteine?
  ord <- do.call(order, args = c(as.list(RES_onoff2[, validvalue_cols]), decreasing = TRUE))
  #cl <- hclust(dist(D_onoff_wide2[, cols], method = "manhattan"), method="complete")


  RES_onoff2_long[, protein_name_column] <- factor(RES_onoff2_long[, protein_name_column],
                                     levels = RES_onoff2[, protein_name_column][ord])

<<<<<<< HEAD
  pl <- ggplot2::ggplot(data = RES_onoff2_long, ggplot2::aes(x = RES_onoff2_long[["group"]], y = RES_onoff2_long[["Gene.names"]], fill = RES_onoff2_long[["value"]])) +  ## TODO: Gene.names
=======
  group <- Gene.names <- value <- NULL # silence notes when checking the package

  pl <- ggplot2::ggplot(data = RES_onoff2_long, ggplot2::aes(x = group, y = Gene.names, fill = value)) +  ## TODO: Gene.names
>>>>>>> cdea481 (clean up to get rid of most warnings and notes when checking)
    ggplot2::geom_tile() +  ggplot2::ylab("Gene name") + ggplot2::xlab("group") + ggplot2::theme_bw()

  #if (onoffGreaterThanEqual < 1 | !is.null(onoffdiff)) {
  pl <- pl + ggplot2::scale_fill_gradient(limits = c(0,max(RES_onoff2_long$value)), low = "white", high = "forestgreen") #
  pl <- pl + ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1.8)),
<<<<<<< HEAD
                            axis.title = ggplot2::element_text(size = ggplot2::rel(1.8)),
                            legend.title = ggplot2::element_text(size= ggplot2::rel(1.8)),
                            legend.text = ggplot2::element_text(size= ggplot2::rel(1.8)))
=======
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1.8)),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(1.8)),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(1.8)))
  pl
>>>>>>> cdea481 (clean up to get rid of most warnings and notes when checking)

  return(pl)

}





