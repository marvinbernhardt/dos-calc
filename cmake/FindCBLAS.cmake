# Find CBLAS
#
# Usage:
#   find_package(CBLAS [REQUIRED])
#
# It sets the following variables:
#   CBLAS_FOUND
#   CBLAS_INCLUDES
#   CBLAS_LIBRARIES

find_path(CBLAS_INCLUDES
    NAMES cblas.h
    PATH_SUFFIXES include inc include/x86_64 include/x64 openblas/include include/blis blis/include blis/include/blis  Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers)
find_library(CBLAS_LIBRARIES
    NAMES cblas blas blis openblas accelerate
    PATH_SUFFIXES lib lib64 lib/x86_64 lib/x64 lib/x86 lib/Win32 lib/import lib64/import)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG
    CBLAS_INCLUDES CBLAS_LIBRARIES)

mark_as_advanced(CBLAS_INCLUDES CBLAS_LIBRARIES)
