#
#  dwave_sapi_FOUND - system has dwave_sapi
#  dwave_sapi_INCLUDE_DIR - the dwave_sapi include directories
#  dwave_sapi_LIBRARY - link these to use dwave_sapi

find_path(dwave_sapi_INCLUDE_DIR NAMES dwave_sapi.h PATHS /usr/local/include)
find_library(dwave_sapi_LIBRARY NAMES libdwave_sapi.dylib PATHS /usr/local/lib)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(dwave_sapi  DEFAULT_MSG
        dwave_sapi_LIBRARY dwave_sapi_INCLUDE_DIR)
