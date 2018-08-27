cmake_minimum_required(VERSION 3.6.2)

project(CBench)

set (CMAKE_CXX_STANDARD 11)

# add libs
set (EXTRA_LIBS "")

# Blosc will always be enabled now
set (CBENCH_ENABLE_BLOSC 1)
#option(CBENCH_ENABLE_BLOSC
#	"Enable BLOSC."
#	ON
#)
option(CBENCH_ENABLE_BIG_CRUNCH
	"Enable BigCrunch."
	ON
)
option(CBENCH_ENABLE_SZ
	"Enable SZ."
	ON
)

option(CBENCH_ENABLE_NYX
	"Enables the reading of NyX HDF5 data"
	ON
)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/metrics)
include_directories(${CMAKE_SOURCE_DIR}/utils)

include_directories(${CMAKE_SOURCE_DIR}/dataLoader)
include_directories(${CMAKE_SOURCE_DIR}/dataLoader/HACC)
if(CBENCH_ENABLE_NYX)
	include_directories(${CMAKE_SOURCE_DIR}/dataLoader/NYX)
endif()

include_directories(${CMAKE_SOURCE_DIR}/thirdparty/genericio)
if(CBENCH_ENABLE_NYX)
    set(HDF5_USE_STATIC_LIBRARIES ON)
	find_package(HDF5 COMPONENTS CXX REQUIRED)
	if(HDF5_FOUND)
		include_directories(${HDF5_INCLUDE_DIRS})
		set (EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARIES} )
		add_definitions (-DCBENCH_HAS_NYX)
	endif()
endif()


# HACC Params
add_definitions(-DGENERICIO_NO_MPI)
add_definitions(-DGENERICIO_NO_COMPRESSION)
add_library(GenericIO STATIC "${CMAKE_SOURCE_DIR}/thirdparty/genericio/GenericIO.cxx")
set_property(TARGET GenericIO PROPERTY POSITION_INDEPENDENT_CODE ON)


# Add Compressors
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)

include_directories(${CMAKE_SOURCE_DIR}/compressors)

if(CBENCH_ENABLE_BLOSC)
    include_directories(${CMAKE_SOURCE_DIR}/compressors/BLOSC)
    if (NOT CBENCH_ENABLE_BIG_CRUNCH)
        add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/c-blosc-1.10.2)
        include_directories(${CMAKE_SOURCE_DIR}/thirdparty/c-blosc-1.10.2/blosc)
        set (EXTRA_LIBS ${EXTRA_LIBS} blosc)
		message(STATUS "--Using external BLOSC--")
    endif()
endif()

if(CBENCH_ENABLE_BIG_CRUNCH)
    include_directories(${CMAKE_SOURCE_DIR}/compressors/BigCrunch)
    add_subdirectory("${CMAKE_SOURCE_DIR}/thirdparty/BigCrunch")
    include_directories(${CMAKE_SOURCE_DIR}/thirdparty/BigCrunch/include)
    add_definitions (-DCBENCH_HAS_BIG_CRUNCH)
    set (EXTRA_LIBS ${EXTRA_LIBS} bigcrunch)
	set (EXTRA_LIBS ${EXTRA_LIBS} blosc)
	message(STATUS "--Using BigCrunchs BLOSC--")
endif()

if(CBENCH_ENABLE_SZ)
    include_directories(${CMAKE_SOURCE_DIR}/compressors/SZ)
    add_subdirectory("${CMAKE_SOURCE_DIR}/thirdparty/SZ")
    include_directories (
        ${CMAKE_SOURCE_DIR}/thirdparty/SZ/sz/include
        ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zlib
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/common
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/compress
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/decompress
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/deprecated
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/dictBuilder
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/dll
            ${CMAKE_SOURCE_DIR}/thirdparty/SZ/zstd/legacy
            )
    add_definitions (-DCBENCH_HAS_SZ)
    set (EXTRA_LIBS ${EXTRA_LIBS} SZ zstd zlib)
endif()


# MPI
find_package(MPI)
if (MPI_FOUND)
	include_directories(${MPI_INCLUDE_PATH})
	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
	SET(CMAKE_C_COMPILER mpicc)
	SET(CMAKE_CXX_COMPILER mpicxx)
else()
	message(FATAL_ERROR "Could not find any MPI installation!!!") 
endif()


# Build executable
add_executable(CBench main.cpp)
target_link_libraries(CBench ${MPI_LIBRARIES} ${EXTRA_LIBS} GenericIO)


# Prevent in source builds since they are evil
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
	message(FATAL_ERROR "In source build is not supported, please restore modified files, 
						delete CMake generated things (CMakeFiles folder and CMakeCache.txt)
						, and build out of source")	
endif()


# Remove output files
add_custom_target(clean-output
	COMMAND rm -f *.log
)
