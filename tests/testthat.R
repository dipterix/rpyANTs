# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(rpyANTs)


message("Checking if rpymat is set up")
if(rpyANTs:::rpymat_is_setup()) {

  message("Checking if ants is available")
  if( rpyANTs:::ants_available() ) {

    message("ANTsPy is set up, run tests...")
    test_check("rpyANTs")

  }
}
