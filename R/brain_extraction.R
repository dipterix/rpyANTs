#' @title Extract brain and strip skull
#' @description
#' Print \code{antspynet$brain_extraction} to see the original documentation.
#' @param x input image or image path
#' @param modality modality type
#' @param verbose whether to print out process to the screen
#' @returns Brain mask image
#' @export
antspynet_brain_extraction <- function(
    x, modality = c("t1", "t1nobrainer", "t1combined", "flair", "t2", "t2star", "bold", "fa", "t1t2infant", "t1infant", "t2infant"),
    verbose=FALSE) {

  antspynet <- load_antspynet()
  image <- as_ANTsImage(x)
  modality <- convert_if_not_python(modality, { match.arg(modality) })
  verbose <- convert_if_not_python(verbose, { isTRUE(as.logical(verbose)) })
  re <- antspynet$brain_extraction(image = image, modality = modality, verbose = verbose)
  return(re)
}
