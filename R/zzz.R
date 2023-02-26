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

  get_ants <- function(force = FALSE) {
    if(!force && inherits(ants, "python.builtin.module")) {
      return( ants )
    }
    if( !rpymat_is_setup() ) {
      return( NULL )
    }
    tryCatch({
      rpymat::ensure_rpymat(verbose = FALSE)
      m <- reticulate::import("ants", convert = FALSE, delay_load = FALSE)
      m <- inject_ants(m)
      class(m) <- c('ants.proxy', class(m))
      ants <<- m
      return( ants )
    }, error = function(e) {
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
    fun = function() {
      load_ants()
    }
  )
  makeActiveBinding("py", fun = load_py, env = pkg)
}

.onUnload <- function(libpath) {

  .ants$clean()
}

#' @title Check if 'ANTs' is available
#' @seealso \code{\link{install_ants}}
#' @returns Logical, whether 'ANTs' is available
#' @export
ants_available <- function() {
  if( !rpymat_is_setup() ) {
    return( FALSE )
  }
  tryCatch({
    rpymat::ensure_rpymat(verbose = FALSE)
    return(reticulate::py_module_available("ants"))
  }, error = function(e) {
    return(FALSE)
  })
}

#' @name ants
#' @title Get 'ANTs' module
#' @param force whether to force reloading \code{ants} module; default is false
#' @usage ants
#' @returns A 'Python' module
#' @export
NULL

#' @rdname ants
#' @export
load_ants <- .ants$get

#' @name py
#' @title Get 'Python' main process environment
#' @usage py
#' @returns The 'Python' main process as a module
#' @export
"py"
NULL

