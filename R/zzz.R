
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
      m <- reticulate::import("ants")
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
#' @export
ants_available <- function() {
  ants <- load_ants()
  !is.null(ants)
}

#' @name ants
#' @title Get 'ANTs' module
#' @param force whether to force reloading \code{ants} module; default is false
#' @usage ants
#' @export
NULL

#' @rdname ants
#' @export
load_ants <- .ants$get

#' @name py
#' @title Get 'Python' built-in main module
#' @usage py
#' @export
"py"
NULL
