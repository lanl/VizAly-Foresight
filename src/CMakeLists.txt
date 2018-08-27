cmake_minimum_required(VERSION 3.6.2)

project(CBench)

set (CMAKE_CXX_STANDARD 11)


include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/metrics)
include_directories(${CMAKE_SOURCE_DIR}/utils)

include_directories(${CMAKE_SOURCE_DIR}/dataLoader)
include_directories(${CMAKE_SOURCE_DIR}/dataLoader/HACC)
# something for nyx here

include_directories(${CMAKE_SOURCE_DIR}/thirdparty/genericio)
# something for nyx here


# HACC Params
add_definitions(-DGENERICIO_NO_MPI)
add_definitions(-DGENERICIO_NO_COMPRESSION)
add_library(GenericIO STATIC "${CMAKE_SOURCE_DIR}/thirdparty/genericio/GenericIO.cxx")
set_property(TARGET GenericIO PROPERTY POSITION_INDEPENDENT_CODE ON)


# add libs
set (EXTRA_LIBS "")

include_directories(${CMAKE_SOURCE_DIR}/compressors)

option(CBENCH_ENABLE_BLOSC "Enable BLOSC." ON)
if(CBENCH_ENABLE_BLOSC)
    include_directories(${CMAKE_SOURCE_DIR}/compressors/BLOSC)

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(BLOSC_LIBRARY ${BLOSC_LIBRARY} CACHE STRING "BLOSC_LIBRARY library")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${BLOSC_LIBRARY})

    #Asks for zfp include files + includes it in the project
    set (BLOSC_INCLUDE_PATH ${BLOSC_INCLUDE_PATH} CACHE STRING "BLOSC include directory")
    include_directories(${BLOSC_INCLUDE_PATH})

    add_definitions (-DCBENCH_HAS_BLOSC)
endif()


option(CBENCH_ENABLE_BIG_CRUNCH "Enable BigCrunch." OFF)
if(CBENCH_ENABLE_BIG_CRUNCH)
    include_directories(${CMAKE_SOURCE_DIR}/compressors/BigCrunch)

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(BIGCRUNCH_LIBRARY ${BIGCRUNCH_LIBRARY} CACHE STRING "BIGCRUNCH_LIBRARY library")
    find_library(BIGCRUNCH_BLOSC_LIBRARY ${BIGCRUNCH_BLOSC_LIBRARY} CACHE STRING "BIGCRUNCH_LIBRARY library")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${BIGCRUNCH_LIBRARY} ${BIGCRUNCH_BLOSC_LIBRARY})

    #Asks for zfp include files + includes it in the project
    set (BIGCRUNCH_INCLUDE_PATH ${BIGCRUNCH_INCLUDE_PATH} CACHE STRING "BIGCRUNCH include directory")
    include_directories(${BIGCRUNCH_INCLUDE_PATH})

    add_definitions (-DCBENCH_HAS_BIG_CRUNCH)
    set (EXTRA_LIBS ${EXTRA_LIBS} ${BIGCRUNCH_LIBRARY} -fopenmp)
endif()


option(CBENCH_ENABLE_SZ "Enable SZ." ON)
if(CBENCH_ENABLE_SZ)
    include_directories(${CMAKE_SOURCE_DIR}/compressors/SZ)

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(SZ_LIBRARY ${SZ_LIBRARY} CACHE STRING "SZ library")
    find_library(ZLIB_LIBRARY ${ZLIB_LIBRARY} CACHE STRING "ZLIB library")
    find_library(ZSTD_LIBRARY ${ZSTD_LIBRARY} CACHE STRING "SZ library")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${SZ_LIBRARY} ${ZLIB_LIBRARY} ${ZSTD_LIBRARY})

    #Asks for zfp include files + includes it in the project
    set (SZ_INCLUDE_PATH ${SZ_INCLUDE_PATH} CACHE STRING "SZ include directory")
    include_directories(${SZ_INCLUDE_PATH})

    add_definitions (-DCBENCH_HAS_SZ)
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
