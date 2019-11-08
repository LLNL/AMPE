#[==[
Adapted from
  * https://github.com/Kitware/VTK/blob/master/CMake/FindNetCDF.cmake
    - Downloaded 10/29/2019
  * https://github.com/lammps/lammps/blob/master/cmake/Modules/FindNetCDF.cmake
    - Downloaded 11/06/2019

Provides the following variables:

  * `NetCDF_FOUND`: Whether NetCDF was found or not.
  * `NetCDF_INCLUDE_DIRS`: Include directories necessary to use NetCDF.
  * `NetCDF_LIBRARIES`: Libraries necessary to use NetCDF.
  * `NetCDF_VERSION`: The version of NetCDF found.
  * `NetCDF::NetCDF`: A target to use with `target_link_libraries`.

You can require certain interfaces to be FOUND by passing a COMPONENTS argument
to find_package containing one or more of

  * CXX
  * CXX_LEGACY
  * F77
  * F90

When interfaces are requested the user has access to interface specific hints:

  * 'NetCDF_${LANG}_INCLUDE_DIR': where to search for interface header files
  * 'NetCDF_${LANG}_LIBRARY': interface library

and the following target(s) will be created

  * 'NetCDF::NetCDF_${LANG}'

#]==]

include(FindPackageHandleStandardArgs)

if(NOT NetCDF_FOUND)
  # Try to find a CMake-built NetCDF.
  message(STATUS "Searching for NetCDF...")
  find_package(netCDF CONFIG QUIET)
  if (netCDF_FOUND)
    # Forward the variables in a consistent way.
    set(NetCDF_FOUND "${netCDF_FOUND}")
    set(NetCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")
    set(NetCDF_LIBRARIES "${netCDF_LIBRARIES}")
    set(NetCDF_LIBRARY_DIRS "${netCDF_LIB_DIR}")
    set(NetCDF_VERSION "${NetCDFVersion}")
    if (NOT TARGET NetCDF::NetCDF)
      add_library(NetCDF::NetCDF INTERFACE IMPORTED)
      set_target_properties(NetCDF::NetCDF PROPERTIES
        INTERFACE_LINK_LIBRARIES "${NetCDF_LIBRARIES}")
    endif ()
  endif ()
endif()

if(NOT NetCDF_FOUND)
  message(STATUS "Searching for PkgConfig...")
  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    pkg_check_modules(_NetCDF QUIET netcdf IMPORTED_TARGET)
    if (_NetCDF_FOUND)
      # Forward the variables in a consistent way.
      set(NetCDF_FOUND "${_NetCDF_FOUND}")
      set(NetCDF_INCLUDE_DIRS "${_NetCDF_INCLUDE_DIRS}")
      set(NetCDF_LIBRARIES "${_NetCDF_LIBRARIES}")
      set(NetCDF_VERSION "${_NetCDF_VERSION}")
      if (NOT TARGET NetCDF::NetCDF)
        add_library(NetCDF::NetCDF INTERFACE IMPORTED)
        set_target_properties(NetCDF::NetCDF PROPERTIES
          INTERFACE_LINK_LIBRARIES "PkgConfig::_NetCDF")
      endif ()
      message(STATUS "NetCDF_LIBRARIES: ${NetCDF_LIBRARIES}")
    endif ()
  endif ()
endif()

if(NOT NetCDF_FOUND)
  message(STATUS "Searching for path to NetCDF include...")
  find_path(NetCDF_INCLUDE_DIR
    NAMES netcdf.h
    DOC "netcdf include directories")
  mark_as_advanced(NetCDF_INCLUDE_DIR)

  find_library(NetCDF_LIBRARY
    NAMES netcdf
    DOC "netcdf library")
  mark_as_advanced(NetCDF_LIBRARY)

  if (NetCDF_INCLUDE_DIR)
    file(STRINGS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
      REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
    string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1"
        _netcdf_version_major "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1"
        _netcdf_version_minor "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1"
        _netcdf_version_patch "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1"
        _netcdf_version_note "${_netcdf_version_lines}")
    set(NetCDF_VERSION
        "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
    unset(_netcdf_version_major)
    unset(_netcdf_version_minor)
    unset(_netcdf_version_patch)
    unset(_netcdf_version_note)
    unset(_netcdf_version_lines)
  endif ()

  find_package_handle_standard_args(NetCDF
    REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR
    VERSION_VAR NetCDF_VERSION)

  if (NetCDF_FOUND)
    set(NetCDF_INCLUDE_DIRS "${NetCDF_INCLUDE_DIR}")
    set(NetCDF_LIBRARIES "${NetCDF_LIBRARY}")
    get_filename_component(NetCDF_LIBRARY_DIRS "${NetCDF_LIBRARY}" DIRECTORY)
    message(STATUS "NetCDF libraries: ${NetCDF_LIBRARIES}")
    if (NOT TARGET NetCDF::NetCDF)
      add_library(NetCDF::NetCDF UNKNOWN IMPORTED)
      set_target_properties(NetCDF::NetCDF PROPERTIES
        IMPORTED_LOCATION "${NetCDF_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIR}")
    endif ()
  endif ()
endif()


# Find requested language components

macro(NetCDF_check_interface lang header libs)
  if(NetCDF_${lang})
    #search starting from user modifiable cache var
    find_path(NetCDF_${lang}_INCLUDE_DIR NAMES ${header}
      HINTS "${NetCDF_INCLUDE_DIRS}"
      HINTS "${NetCDF_${lang}_ROOT}/include"
      ${USE_DEFAULT_PATHS})

    find_library(NetCDF_${lang}_LIBRARY NAMES ${libs}
      HINTS "${NetCDF_lib_dirs}"
      HINTS "${NetCDF_${lang}_ROOT}/lib"
      ${USE_DEFAULT_PATHS})

    mark_as_advanced(NetCDF_${lang}_INCLUDE_DIR NetCDF_${lang}_LIBRARY)

    #export to internal varS that rest of project can use directly
    set(NetCDF_${lang}_LIBRARIES ${NetCDF_${lang}_LIBRARY})
    set(NetCDF_${lang}_INCLUDE_DIRS ${NetCDF_${lang}_INCLUDE_DIR})

    find_package_handle_standard_args(NetCDF
      REQUIRED_VARS NetCDF_${lang}_LIBRARY NetCDF_${lang}_INCLUDE_DIR)

    list(APPEND NetCDF_LIBRARIES ${NetCDF_${lang}_LIBRARY})
    list(APPEND NetCDF_INCLUDE_DIRS ${NetCDF_${lang}_INCLUDE_DIR})

    if(NOT TARGET NetCDF::NetCDF_${lang})
      add_library(NetCDF::NetCDF_${lang} INTERFACE IMPORTED)
      set_property(TARGET NetCDF::NetCDF_${lang} PROPERTY
          INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_${lang}_INCLUDE_DIR}")
      set_property(TARGET NetCDF::NetCDF_${lang} APPEND PROPERTY
          INTERFACE_LINK_LIBRARIES "${NetCDF_${lang}_LIBRARY}" NetCDF::NetCDF)
    endif()
  endif()
endmacro()

list(FIND NetCDF_FIND_COMPONENTS "CXX" _nextcomp)
if(_nextcomp GREATER -1)
  set(NetCDF_CXX 1)
endif()
list(FIND NetCDF_FIND_COMPONENTS "CXX_LEGACY" _nextcomp)
if(_nextcomp GREATER -1)
  set(NetCDF_CXX_LEGACY 1)
endif()
list(FIND NetCDF_FIND_COMPONENTS "F77" _nextcomp)
if(_nextcomp GREATER -1)
  set(NetCDF_F77 1)
endif()
list(FIND NetCDF_FIND_COMPONENTS "F90" _nextcomp)
if(_nextcomp GREATER -1)
  set(NetCDF_F90 1)
endif()
NetCDF_check_interface(CXX netcdf netcdf_c++4)
NetCDF_check_interface(CXX_LEGACY netcdfcpp.h netcdf_c++)
NetCDF_check_interface(F77 netcdf.inc netcdff)
NetCDF_check_interface(F90 netcdf.mod netcdff)

list(REMOVE_DUPLICATES NetCDF_INCLUDE_DIRS)
list(REMOVE_DUPLICATES NetCDF_LIBRARIES)
