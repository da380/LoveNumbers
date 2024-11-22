#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "LoveNumbers::LoveNumbers" for configuration ""
set_property(TARGET LoveNumbers::LoveNumbers APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(LoveNumbers::LoveNumbers PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libLoveNumbers.so"
  IMPORTED_SONAME_NOCONFIG "libLoveNumbers.so"
  )

list(APPEND _cmake_import_check_targets LoveNumbers::LoveNumbers )
list(APPEND _cmake_import_check_files_for_LoveNumbers::LoveNumbers "${_IMPORT_PREFIX}/lib/libLoveNumbers.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
