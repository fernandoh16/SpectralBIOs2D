# NTL needs GMP 3.1.1 or newer, this script will fail if an old version is
# detected

find_package( GMP REQUIRED )

if( NOT GMP_FOUND )

  message( ERROR "NTL requires GMP" )
  set( NTL_FOUND FALSE )

else( NOT GMP_FOUND )

  if( _IS_GMP_VERSION_TOO_LOW )

    message( ERROR, "NTL needs GMP>=3.1.1. Your GMP version is ${CGAL_GMP_VERSION}." )

  else( _IS_GMP_VERSION_TOO_LOW )

    find_path(NTL_INCLUDE_DIR
              NAMES NTL/GF2X.h
              PATHS ENV NTL_INC_DIR
                    ENV NTL_DIR
              PATH_SUFFIXES include
              DOC "The directory containing the NTL header files"
             )

    find_library(NTL_LIBRARY
                 NAMES ntl
                 PATHS ENV NTL_LIB_DIR
                       ENV NTL_DIR
                       ENV LD_LIBRARY_PATH
                 PATH_SUFFIXES lib
                 DOC "Path to the NTL library"
                 /usr/local/lib/libntl.a
                )

    # TODO if NTL_INC_DIR is given you should not search in default path

#    find_library(NTL_LIBRARY
#                 NAMES ntl
#                 PATHS ENV NTL_LIB_DIR
#                 DOC "Path to the NTL library"
#                )

    message( STATUS "NTL_INCLUDE_DIR = '${NTL_INCLUDE_DIR}'" )
    message( STATUS "NTL_LIBRARY = '${NTL_LIBRARY}'" )

    if ( NTL_INCLUDE_DIR AND NTL_LIBRARY ) 
      
       #check version

       set( NTL_VERSION_H "${NTL_INCLUDE_DIR}/NTL/version.h" )

       if ( EXISTS ${NTL_VERSION_H} )

         file( READ "${NTL_VERSION_H}" NTL_VERSION_H_CONTENTS )
  
         string( REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" CGAL_NTL_VERSION "${NTL_VERSION_H_CONTENTS}" )
  
         message( STATUS "Found NTL version: ${CGAL_NTL_VERSION}" )

         #IS_VERSION_GREATER( "${CGAL_NTL_VERSION}" "5.0.0" _IS_NTL_VERSION_GREATER )

         if ( ${CGAL_NTL_VERSION} VERSION_GREATER "5.0.0" )
           set( NTL_FOUND TRUE )
         else ()
           set( NTL_FOUND FALSE )
         endif ()


       endif (EXISTS ${NTL_VERSION_H} )

    endif ( NTL_INCLUDE_DIR AND NTL_LIBRARY ) 

    if ( NTL_FOUND )

      set ( NTL_INCLUDE_DIRS ${NTL_INCLUDE_DIR} )
      set ( NTL_LIBRARIES ${NTL_LIBRARY} )

      get_filename_component(NTL_LIBRARIES_DIR ${NTL_LIBRARIES} PATH CACHE )

      mark_as_advanced( NTL_INCLUDE_DIR NTL_LIBRARY )

    endif( NTL_FOUND )


  endif( _IS_GMP_VERSION_TOO_LOW )

endif( NOT GMP_FOUND )

if ( NTL_FOUND )
else ( NTL_FOUND )
  if ( NTL_FIND_REQUIRED )
    message( FATAL_ERROR "Could not find NTL" )
  endif ( NTL_FIND_REQUIRED )
endif ( NTL_FOUND )
