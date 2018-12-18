cmake_minimum_required(VERSION 3.6.2)

project(CBench)

set (CMAKE_CXX_STANDARD 11)

# add libs
set (EXTRA_LIBS "")

option(CBENCH_ENABLE_NYX_LOADER
	"Enables the reading of NyX HDF5 data"
	OFF
)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/metrics)
include_directories(${CMAKE_SOURCE_DIR}/utils)

include_directories(${CMAKE_SOURCE_DIR}/dataLoader)
include_directories(${CMAKE_SOURCE_DIR}/dataLoader/Binary)
include_directories(${CMAKE_SOURCE_DIR}/dataLoader/HACC)
if(CBENCH_ENABLE_NYX_LOADER)
	include_directories(${CMAKE_SOURCE_DIR}/dataLoader/NYX)
endif()

# HACC Loader
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/genericio)
add_definitions(-DGENERICIO_NO_COMPRESSION)
add_library(GenericIO STATIC "${CMAKE_SOURCE_DIR}/thirdparty/genericio/GenericIO.cxx")
set_property(TARGET GenericIO PROPERTY POSITION_INDEPENDENT_CODE ON)

# NYX Loader
if(CBENCH_ENABLE_NYX_LOADER)
    set(HDF5_USE_STATIC_LIBRARIES ON)
	find_package(HDF5 COMPONENTS CXX REQUIRED)
	if(HDF5_FOUND)
		include_directories(${HDF5_INCLUDE_DIRS})
		set (EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARIES} )
		add_definitions (-DCBENCH_HAS_NYX)
	endif()
endif()
    

# Creates header file, compressorIncludes.h, with all available compressors
file(GLOB_RECURSE header_files "${CMAKE_SOURCE_DIR}/compressors/*/*.hpp")
foreach(f ${header_files})
    file(RELATIVE_PATH fr "${CMAKE_CURRENT_SOURCE_DIR}/compressors/" ${f})
    list(APPEND incls "#include \"${fr}\"\n")
endforeach()
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/compressors/compressorIncludes.h ${incls})


# Add external compressors
include_directories(${CMAKE_SOURCE_DIR}/compressors)
file(GLOB subdir "${CMAKE_SOURCE_DIR}/compressors/*")
foreach (dir ${subdir})
    # Skip files ending in .h
    if ("${dir}" MATCHES "\\.h$")
        continue()
    endif()

    # Skip files ending in .hpp
    if ("${dir}" MATCHES "\\.hpp$")
        continue()
    endif()
    
    # Skip files ending in .md
    if ("${dir}" MATCHES "\\.md$")
        continue()
    endif()


    include_directories(${dir})

    if (NOT EXISTS ${dir}/CMakeLists.txt)
        continue()
    else()
        # Only process folders with a CMakelists.txt
        set_property(GLOBAL PROPERTY includeDirectories "")
        set_property(GLOBAL PROPERTY linkLibraries "")
        set_property(GLOBAL PROPERTY definitions "")
        
        add_subdirectory(${dir})

        #include directories to add
        get_property(incDir GLOBAL PROPERTY includeDirectories)
        include_directories(${incDir})

        #link libraries to add directories to add
        get_property(linkLibs GLOBAL PROPERTY linkLibraries)
        set(EXTRA_LIBS ${EXTRA_LIBS}  ${linkLibs})

        # add defintions
        get_property(def GLOBAL PROPERTY definitions)
        add_definitions (${def})
    endif()
endforeach()



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
