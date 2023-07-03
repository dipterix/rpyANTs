#' @importFrom reticulate tuple
#' @importFrom reticulate py_help
#' @importFrom reticulate np_array
#' @importFrom reticulate import
#' @importFrom reticulate py_bool
#' @importFrom reticulate py_call
#' @importFrom reticulate py_del_attr
#' @importFrom reticulate py_dict
#' @importFrom reticulate py_del_item
#' @importFrom reticulate py_eval
#' @importFrom reticulate py_get_attr
#' @importFrom reticulate py_get_item
#' @importFrom reticulate py_len
#' @importFrom reticulate py_none
#' @importFrom reticulate py_set_attr
#' @importFrom reticulate py_set_item
#' @importFrom reticulate py_str
#' @importFrom reticulate py_to_r
#' @importFrom reticulate py_to_r_wrapper
#' @importFrom reticulate r_to_py
#' @importFrom reticulate py_run_string
#' @importFrom rpymat run_script
#' @importFrom rpymat repl_python
NULL

AFFINE_TRANSFORM_TYPES <- c(
  'AffineTransform', 'CenteredAffineTransform',
  'Euler2DTransform', 'Euler3DTransform', 'Rigid3DTransform',
  'Rigid2DTransform', 'QuaternionRigidTransform',
  'Similarity2DTransform', 'CenteredSimilarity2DTransform',
  'Similarity3DTransform', 'CenteredRigid2DTransform',
  'CenteredEuler3DTransform', "ScaleTransform",
  "ScaleVersor3DTransform", "ScaleSkewVersor3DTransform"
)

validate_python <- function(verbose = TRUE) {
  verb <- function(expr) {
    if(verbose) {
      force( expr )
    }
  }
  verb(message("Initializing python environment: "))

  rpymat::ensure_rpymat(verbose = verbose)

  verb(message("Trying to get installed packages..."))
  tbl <- rpymat::list_pkgs()
  pkgs <- tbl$package
  pkgs <- pkgs[grepl("^[a-zA-Z0-9]", pkgs)]

  verb(cat("Installed packages:", paste(pkgs, collapse = ", "), "\n"))

  # Check environment
  verb(message("Trying to validate packages..."))

  package_missing <- NULL
  for(package in c("numpy", "h5py", "cython", "pandas", "ants")) {
    tryCatch({
      verb({ cat(sprintf("%s: ...", package)) })
      module <- reticulate::import(package)
      verb({ cat("\b\b\b", module$`__version__`, "\n", sep = "") })
    }, error = function(e) {
      verb({ cat("\b\b\bN/A\n", sep = "") })
      package_missing <<- c(package_missing, package)
    })
  }

  return(invisible(package_missing))
}

#' @title Install 'ANTs' via 'ANTsPy'
#' @param python_ver 'Python' version, see \code{\link[rpymat]{configure_conda}};
#' default is \code{"3.9"} since \code{'ANTsPy'} is compiled for all
#' @param verbose whether to print the installation messages
#' @returns This function returns nothing.
#' @export
install_ants <- function(python_ver = "3.9", verbose = TRUE) {
  # Install conda and create a conda environment
  if(!dir.exists(rpymat::env_path())) {
    standalone <- !file.exists(rpymat::conda_bin())
    rpymat::configure_conda(python_ver = python_ver, force = TRUE, standalone = standalone)
  }
  rpymat::ensure_rpymat(verbose = verbose)
  installed_pkgs_tbl <- rpymat::list_pkgs()

  # install necessary libraries
  pkgs <- c("h5py", "numpy", "scipy", "pandas", "cython")
  if(!all(pkgs %in% installed_pkgs_tbl$package)) {
    rpymat::add_packages(pkgs)
  }

  # install antspyx family
  ants_packages <- c("antspyx", "antspynet")
  ants_packages <- ants_packages[!ants_packages %in% installed_pkgs_tbl$package]
  if(length(ants_packages)) {
    rpymat::add_packages(packages = ants_packages, pip = TRUE)
  }

  validate_python(verbose = verbose)
}

rpymat_is_setup <- function() {
  return( dir.exists(rpymat::env_path()) )
}

