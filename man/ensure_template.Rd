% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/templates.R
\name{ensure_template}
\alias{ensure_template}
\title{Ensure the template directory is downloaded}
\usage{
ensure_template(name = BUILTIN_TEMPLATES)
}
\arguments{
\item{name}{name of the template, commonly known as \code{'MNI152'}
templates; choices are \code{"mni_icbm152_nlin_asym_09a"},
\code{"mni_icbm152_nlin_asym_09b"}, and \code{"mni_icbm152_nlin_asym_09c"}.}
}
\value{
The downloaded template path
}
\description{
Ensure the template directory is downloaded
}
\examples{

# Do not run for testing as this will download the template
if(FALSE) {

# Default is `mni_icbm152_nlin_asym_09a`
ensure_template()
ensure_template("mni_icbm152_nlin_asym_09a")

# Using MNI152b
ensure_template("mni_icbm152_nlin_asym_09b")

}


}
