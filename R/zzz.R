#
# ex <- c("n3_bias_field_correction2", "n4_bias_field_correction")
# for(nm in names(ants)) {
#   if(nm %in% ex) { next }
#   message(nm)
#   reticulate::py_help_handler("completion", nm, source = "ants")
# }

inject_ants <- function(ants) {

  if(
    !inherits(ants, "python.builtin.module") ||
    !identical(get_os(), "windows")
  ) { return(ants) }

  ants$utils$rpyANTsInjected <- FALSE

  RpyANTsInjection <- reticulate::PyClass(
    "RpyANTsInjection",
    list(
      requested_tempfiles = NULL,
      `__init__` = function(self) {
        self$requested_tempfiles <- py_list()
        NULL
      },
      temp_nii_file = function(self) {
        f <- tempfile(fileext = ".nii.gz")
        f <- normalizePath(f, mustWork = FALSE, winslash = "/")
        self$requested_tempfiles[[length(self$requested_tempfiles) + 1]] <- f
        return( f )
      }
    )
  )

  ants$utils$rpyANTsInjection <- RpyANTsInjection()

  patch <- system.file(package = "rpyANTs", "patches", "win_patch.py")
  tdir <- normalizePath(tempdir(check = TRUE))
  file.copy(patch, file.path(tdir, "rpyANTs_win_patch.py"))

  reticulate::import_from_path(
    path = tdir,
    module = "rpyANTs_win_patch",
    convert = FALSE,
    delay_load = FALSE
  )
  # ants$utils$win_patch <- win_patch

  ants

}

.ants <- local({

  ants <- NULL

  get_ants <- function(force = FALSE, error_if_missing = TRUE) {
    if(!force && inherits(ants, "python.builtin.module")) {
      return( ants )
    }
    if( !rpymat_is_setup() ) {
      if( error_if_missing ) {
        stop("Please configure ANTsPy environment first. Run the following command:\n  rpyANTs::install_ants()\nIf you would like a specific Python version (e.g. python=3.9), \n  rpyANTs::install_ants(python_ver = '3.9')")
      }
      return( NULL )
    }
    tryCatch({
      rpymat::ensure_rpymat(verbose = FALSE)
      m <- reticulate::import("ants", convert = FALSE, delay_load = FALSE)
      # m <- inject_ants(m)
      class(m) <- c('ants.proxy', class(m))
      ants <<- m
      return( ants )
    }, error = function(e) {
      if( error_if_missing ) {
        stop(e)
      }
      return(NULL)
    })
  }

  clean_ants <- function() {
    ants <<- NULL
  }

  list(
    get = get_ants,
    clean = clean_ants
  )
})

.rpyants <- local({

  rpyants_module <- NULL

  get_rpyants <- function(force = FALSE, error_if_missing = TRUE) {
    if(!force && inherits(rpyants_module, "python.builtin.module")) {
      return( rpyants_module )
    }
    if( !rpymat_is_setup() ) {
      if( error_if_missing ) {
        stop("Please configure ANTsPy environment first. Run the following command:\n  rpyANTs::install_ants()\nIf you would like a specific Python version (e.g. python=3.9), \n  rpyANTs::install_ants(python_ver = '3.9')")
      }
      return( NULL )
    }
    tryCatch({
      rpymat::ensure_rpymat(verbose = FALSE)
      m <- reticulate::import_from_path("rpyants", path = system.file("rpyants", package = "rpyANTs"), convert = FALSE, delay_load = FALSE)
      # m <- reticulate::import("rpyants", convert = FALSE, delay_load = FALSE)
      # m <- inject_ants(m)
      class(m) <- c('ants.proxy', class(m))
      rpyants_module <<- m
      return( rpyants_module )
    }, error = function(e) {
      if( error_if_missing ) {
        stop(e)
      }
      return(NULL)
    })
  }

  clean_rpyants <- function() {
    rpyants_module <<- NULL
  }

  list(
    get = get_rpyants,
    clean = clean_rpyants
  )
})


.antspynet <- local({

  antspynet <- NULL

  get_antspynet <- function(force = FALSE, error_if_missing = TRUE) {
    if(!force && inherits(antspynet, "python.builtin.module")) {
      return( antspynet )
    }
    if( !rpymat_is_setup() ) {
      if( error_if_missing ) {
        stop("Please configure ANTsPy environment first. Run the following command:\n  rpyANTs::install_ants()\nIf you would like a specific Python version (e.g. python=3.9), \n  rpyANTs::install_ants(python_ver = '3.9')")
      }
      return( NULL )
    }
    tryCatch({
      rpyants <- load_rpyants()
      if(is.null(rpyants) || !"try_import_antspynet" %in% names(rpyants$utils$paths)) {
        rpymat::ensure_rpymat(verbose = FALSE)
        m <- reticulate::import("antspynet", convert = FALSE, delay_load = FALSE)
        # set cache directory
        cache_dir <- file_path(R_user_dir(package = "rpyANTs", which = "data"), "keras", "ANTsXNet")
        cache_dir <- dir_create2(cache_dir)
        m$set_antsxnet_cache_directory(cache_dir)
      } else {
        m <- rpyants$utils$paths$try_import_antspynet()
      }
      class(m) <- c('ants.proxy', class(m))
      antspynet <<- m
      return( antspynet )
    }, error = function(e) {
      if( error_if_missing ) {
        stop(e)
      }
      return(NULL)
    })
  }

  clean_antspynet <- function() {
    antspynet <<- NULL
  }

  list(
    get = get_antspynet,
    clean = clean_antspynet
  )
})


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
  makeActiveBinding("py", fun = load_py, env = pkg)
  makeActiveBinding(
    "ants", env = pkg,
    fun = function() {
      load_ants(error_if_missing = FALSE)
    }
  )
  # makeActiveBinding(
  #   "antspynet", env = pkg,
  #   fun = function() {
  #     load_antspynet(error_if_missing = FALSE)
  #   }
  # )
}

.onUnload <- function(libpath) {
  .rpyants$clean()
  .antspynet$clean()
  .ants$clean()
}

#' @title Check if 'ANTs' is available
#' @seealso \code{\link{install_ants}}
#' @param module either \code{'ants'} or \code{'antspynet'}; default is
#' \code{'ants'}
#' @returns Logical, whether \code{'ANTs'} or \code{'ANTsPyNet'} is available
#' @export
ants_available <- function(module = c("ants", "antspynet")) {

  module <- match.arg(module)

  if( !rpymat_is_setup() ) {
    return( FALSE )
  }
  tryCatch({
    rpymat::ensure_rpymat(verbose = FALSE)
    return( reticulate::py_module_available(module) )
  }, error = function(e) {
    return( FALSE )
  })
}

#' @name ants
#' @title Get 'ANTsPy' module
#' @param force whether to force reloading \code{ants} module; default is false
#' @param error_if_missing whether to raise errors when the module is unable to
#' load; default is true.
#' @usage ants
#' @returns A 'Python' module if successfully loaded. If \code{error_if_missing}
#' is set to false and module is unable to load, return \code{NULL}
#' @seealso \code{\link{antspynet}}
#' @export
NULL

#' @rdname ants
#' @export
load_ants <- .ants$get

load_rpyants <- .rpyants$get


#' @name antspynet
#' @title Get \code{'ANTsPyNet'} module
#' @param force whether to force reloading \code{antspynet} module; default is false
#' @param error_if_missing whether to raise errors when the module is unable to
#' load; default is true.
#' @returns A 'Python' module if successfully loaded. If \code{error_if_missing}
#' is set to false and module is unable to load, return \code{NULL}
#' @seealso \code{\link{ants}}
#' @export
load_antspynet <- .antspynet$get

#' @name py
#' @title Get 'Python' main process environment
#' @usage py
#' @returns The 'Python' main process as a module
#' @export
"py"
NULL


