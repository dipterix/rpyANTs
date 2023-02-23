

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

#' @export
install_ants <- function(python_ver = "auto", verbose = TRUE) {
  # Install conda and create a conda environment
  if(!dir.exists(rpymat::env_path())) {
    rpymat::configure_conda(python_ver = python_ver, force = TRUE)
  }
  rpymat::ensure_rpymat(verbose = verbose)
  installed_pkgs_tbl <- rpymat::list_pkgs()

  # install necessary libraries
  pkgs <- c("h5py", "numpy", "scipy", "pandas", "cython")
  if(!all(pkgs %in% installed_pkgs_tbl$package)) {
    rpymat::add_packages(pkgs)
  }

  # install antspyx family
  if(!"antspyx" %in% installed_pkgs_tbl$package) {
    rpymat::add_packages(packages = "antspyx", pip = TRUE)
  }

  validate_python(verbose = verbose)
}

#' @name ants
#' @title Get 'ANTs' module
#' @export
NULL

#' @name py
#' @title Get 'Python' object
#' @export
NULL

rpymat_is_setup <- function() {
  return( dir.exists(rpymat::env_path()) )
}

#' @rdname ants
#' @export
load_ants <- local({

  ants <- NULL

  function(force = FALSE) {
    if(!force && inherits(ants, "python.builtin.module")) {
      return( ants )
    }
    if( !rpymat_is_setup() ) {
      return( NULL )
    }
    tryCatch({
      rpymat::ensure_rpymat(verbose = FALSE)
      ants <<- reticulate::import("ants")
      return( ants )
    }, error = function(e) {
      return(NULL)
    })
  }
})

#' @rdname ants
#' @export
ants_available <- function() {
  ants <- load_ants()
  !is.null(ants)
}

load_py <- local({

  main <- NULL

  function() {
    if (!is.null(main)) { return(main) }

    if( !rpymat_is_setup() ) {
      return( NULL )
    }

    py <- tryCatch({
      reticulate <- asNamespace("reticulate")
      if(isTRUE(reticulate$is_python_initialized())) {
        py <- reticulate::import_main(convert = TRUE)
      } else {
        py <- NULL
      }
      py
    }, error = function(e) {
      reticulate::py
    })

    if(!is.null(py)) {
      main <<- py
    }
    main
  }
})

.onLoad <- function(libname, pkgname) {
  pkg <- getNamespace(pkgname)
  makeActiveBinding(
    "ants", env = pkg,
    fun = function() { load_ants() }
  )
  makeActiveBinding("py", fun = load_py, env = pkg)
}


#' @export
tuple <- reticulate::tuple

#' @export
py_help <- reticulate::py_help

#' @export
np_array <- reticulate::np_array

#' @export
import <- reticulate::import
