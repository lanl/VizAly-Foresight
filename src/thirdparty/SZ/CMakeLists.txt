cmake_minimum_required (VERSION 3.1.3)
project (sz)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm -fopenmp")

set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

find_package (FFTW)
if (FFTW_FOUND)
  include_directories (${FFTW_INCLUDES})
  add_definitions (-DHAVE_FFTW3)
endif ()

find_package (NETCDF)
if (NETCDF_FOUND)
  include_directories (${NETCDF_INCLUDE_DIRS})
endif ()

include_directories (
	${CMAKE_CURRENT_SOURCE_DIR}/sz/include 
	${CMAKE_CURRENT_SOURCE_DIR}/sz/include 
	${CMAKE_CURRENT_SOURCE_DIR}/zlib 
	${CMAKE_CURRENT_SOURCE_DIR}/zstd 
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/common 
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/compress
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/decompress
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/deprecated
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/dictBuilder
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/dll
	${CMAKE_CURRENT_SOURCE_DIR}/zstd/legacy
	)
# include_directories (${ZLIB_INCLUDE_DIRS})

add_subdirectory (sz)
add_subdirectory (zlib)
add_subdirectory (zstd)
add_subdirectory (example)
