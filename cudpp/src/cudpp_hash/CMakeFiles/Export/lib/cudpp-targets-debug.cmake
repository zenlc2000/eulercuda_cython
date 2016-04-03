#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cudpp" for configuration "Debug"
set_property(TARGET cudpp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(cudpp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "/Developer/NVIDIA/CUDA-7.5/lib/libcudart_static.a;-Wl,-rpath,/usr/local/cuda/lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libcudpp64d.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS cudpp )
list(APPEND _IMPORT_CHECK_FILES_FOR_cudpp "${_IMPORT_PREFIX}/lib/libcudpp64d.a" )

# Import target "cudpp_hash" for configuration "Debug"
set_property(TARGET cudpp_hash APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(cudpp_hash PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "/Developer/NVIDIA/CUDA-7.5/lib/libcudart_static.a;-Wl,-rpath,/usr/local/cuda/lib;cudpp"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libcudpp_hash64d.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS cudpp_hash )
list(APPEND _IMPORT_CHECK_FILES_FOR_cudpp_hash "${_IMPORT_PREFIX}/lib/libcudpp_hash64d.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
