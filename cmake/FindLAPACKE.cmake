# Find LAPACKE
#
# Usage:
#   find_package(LAPACKE [REQUIRED])
#
# It sets the following variables:
#   LAPACKE_FOUND
#   LAPACKE_INCLUDES
#   LAPACKE_LIBRARIES

find_path(LAPACKE_INCLUDES
    NAMES lapacke.h
    PATH_SUFFIXES lapack lapacke)
find_library(LAPACKE_LIBRARIES
    NAMES lapacke reflapacke openblas
    PATH_SUFFIXES lapack lapacke)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE DEFAULT_MSG
    LAPACKE_INCLUDES LAPACKE_LIBRARIES)

mark_as_advanced(LAPACKE_INCLUDES LAPACKE_LIBRARIES)
