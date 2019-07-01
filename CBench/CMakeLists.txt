cmake_minimum_required(VERSION 3.6.2)
cmake_policy(SET CMP0048 NEW)

# version actually set for Cori
project(CBench VERSION 1.0.0 LANGUAGES CXX C)

# Prevent in source builds since they are evil
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In source build is not supported, please restore modified files,
                        delete CMake generated things (CMakeFiles folder and CMakeCache.txt)
                        , and build out of source")
endif()

# Set the main executable
add_executable(CBench main.cpp)

set(DEFAULT_BUILD_TYPE "Debug")

# Set build type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()


if(WIN32)
    target_compile_definitions(CBench PRIVATE -DNOMINMAX)
endif()

# Set build options
option(CBENCH_ENABLE_NYX_LOADER "Enables NYX HDF5 data loader" OFF)
option(CBENCH_ENABLE_BIG_CRUNCH "Enable BigCrunch." OFF)
option(CBENCH_ENABLE_BLOSC      "Enable BLOSC." OFF)
option(CBENCH_ENABLE_FPZIP      "Enable FPZIP." OFF)
option(CBENCH_ENABLE_ISABELA    "Enable ISABELA." OFF)
option(CBENCH_ENABLE_LOSSY_WAVE "Enable LossyWave." OFF)
option(CBENCH_ENABLE_MGARD      "Enable MGARD." OFF)
option(CBENCH_ENABLE_SZ         "Enable SZ." OFF)
option(CBENCH_ENABLE_ZFP        "Enable ZFP." OFF)

# include dirs
target_compile_features(CBench PUBLIC cxx_std_14)
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}")
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/metrics")
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/utils")
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/dataLoader")
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/dataLoader/Binary")
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/dataLoader/HACC")


# HACC Loader
target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/thirdparty/genericio")
target_compile_definitions(CBench PRIVATE -DGENERICIO_NO_COMPRESSION)

add_library(GenericIO STATIC)
target_sources(GenericIO PRIVATE "${CMAKE_SOURCE_DIR}/thirdparty/genericio/GenericIO.cxx")
set_target_properties(GenericIO PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(GenericIO PRIVATE -DGENERICIO_NO_COMPRESSION)
target_link_libraries(CBench PRIVATE GenericIO)


# NYX loader
if (CBENCH_ENABLE_NYX_LOADER)
    # NYX Loader
    set(HDF5_USE_STATIC_LIBRARIES ON)

    # assume that CBENCH_PLATFORM is set by the user
    if (${CBENCH_PLATFORM} MATCHES "Cori")
        set(HDF5_CORI "/usr/common/software/hdf5-parallel/1.10.1/gnu")
        find_package(HDF5 REQUIRED PATHS ${HDF5_CORI} NO_DEFAULT_PATH COMPONENTS C)
        target_link_libraries(CBench PRIVATE dl m z hdf5 hdf5_hl)
    else ()
        find_package(HDF5 REQUIRED COMPONENTS C)
        target_link_libraries(CBench PRIVATE ${HDF5_LIBRARIES})
    endif()
    target_include_directories(CBench PRIVATE "${HDF5_INCLUDE_DIRS}")
    target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_NYX)
    target_include_directories(CBench PRIVATE "${CMAKE_SOURCE_DIR}/dataLoader/NYX")
endif()


# Creates header file, compressorIncludes.h, with all available compressors
file(GLOB_RECURSE header_files "${CMAKE_SOURCE_DIR}/compressors/*/*.hpp")
foreach(f ${header_files})
    file(RELATIVE_PATH fr "${CMAKE_CURRENT_SOURCE_DIR}/compressors/" ${f})
    list(APPEND incls "#include \"${fr}\"\n")
endforeach()
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/compressors/compressorIncludes.h ${incls})


# Add external compressors
target_include_directories(CBench PRIVATE ${CMAKE_SOURCE_DIR}/compressors)
file(GLOB SUBDIR "${CMAKE_SOURCE_DIR}/compressors/*")
foreach (COMPRESSOR ${SUBDIR})
    # Only process folders with a CMakelists.txt
    if ((IS_DIRECTORY ${COMPRESSOR}) AND (EXISTS ${COMPRESSOR}/CMakeLists.txt))
        add_subdirectory(${COMPRESSOR})
        get_target_property(CURRENT_LIB CBench LIB)
        if (CURRENT_LIB)
          target_link_libraries(CBench PRIVATE ${CURRENT_LIB})
        endif()
    endif ()
endforeach()


# MPI
find_package(MPI REQUIRED)

# For supporting CMake < 3.9:
if(NOT TARGET MPI::MPI_CXX)
    add_library(MPI::MPI_CXX IMPORTED INTERFACE)
    set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_COMPILE_OPTIONS "${MPI_COMPILE_FLAGS}")
    set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_INCLUDE_PATH}")
    set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES "${MPI_LINK_FLAGS}")
endif()
target_link_libraries(CBench PRIVATE MPI::MPI_CXX)
target_link_libraries(GenericIO PRIVATE MPI::MPI_CXX)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    target_link_libraries(CBench PRIVATE OpenMP::OpenMP_CXX)
endif()


# Remove output files
add_custom_target(clean-output COMMAND rm -f *.log)
