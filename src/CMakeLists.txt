cmake_minimum_required(VERSION 3.6.2)

project(CBench)

set (CMAKE_CXX_STANDARD 11)


include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/metrics)
include_directories(${CMAKE_SOURCE_DIR}/utils)

include_directories(${CMAKE_SOURCE_DIR}/compressors)
include_directories(${CMAKE_SOURCE_DIR}/compressors/BLOSC)
include_directories(${CMAKE_SOURCE_DIR}/compressors/BigCrunch)
include_directories(${CMAKE_SOURCE_DIR}/compressors/SZ)



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


# Add Compressors
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)


add_subdirectory("${CMAKE_SOURCE_DIR}/thirdparty/BigCrunch")
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/BigCrunch/include)

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

# add libs
set (EXTRA_LIBS "")


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
target_link_libraries(CBench ${MPI_LIBRARIES} ${EXTRA_LIBS} GenericIO bigcrunch SZ zstd zlib)


# Prevent in source builds since they are evil
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
	message(FATAL_ERROR "In source build is not supported, please restore modified files, 
						delete CMake generated things (CMakeFiles folder and CMakeCache.txt)
						, and build out of source")	
endif()


# Remove output files
add_custom_target(clean-output
	COMMAND rm -f *.vtu
	COMMAND rm -f *.pvtu
	COMMAND rm -f *.vts
	COMMAND rm -f *.pvts
)
