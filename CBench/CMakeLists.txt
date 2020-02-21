cmake_minimum_required(VERSION 3.6.2)

project(CBench)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel.")

# add libs
set (EXTRA_LIBS "")


# inclue dirs
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/metrics)
include_directories(${CMAKE_SOURCE_DIR}/utils)
include_directories(${CMAKE_SOURCE_DIR}/dataLoader)



# HACC Loader
include_directories(${CMAKE_SOURCE_DIR}/dataLoader/HACC)


# VPIC GDA Loader
option(CBENCH_LOADER_GDA "Enables the reading of VPIC GDA data" OFF)
if (CBENCH_LOADER_GDA)
	add_definitions(-DCBENCH_HAS_GDA)
	include_directories(${CMAKE_SOURCE_DIR}/dataLoader/VPIC_GDA)
endif()


# AMRex Plt File Loader
option(CBENCH_LOADER_AMReX "Enables the reading of AMRex Plt Files" OFF)
if (CBENCH_LOADER_AMReX)
	add_definitions(-DCBENCH_HAS_AMReX)
	set (AMReX_PATH ${AMReX_PATH} CACHE PATH "AMReX directory")
	include_directories(${AMReX_PATH}/include)
	link_directories(${AMReX_PATH}/lib)
	set (EXTRA_LIBS ${EXTRA_LIBS} libamrex.so)
endif()


# VTK Loader
option(CBENCH_LOADER_VTK "Enables the reading of VTK data" OFF)
if (CBENCH_LOADER_VTK)
	find_package(VTK REQUIRED)
	if (VTK_FOUND)
		add_definitions(-DCBENCH_HAS_VTK)
		include_directories(${CMAKE_SOURCE_DIR}/dataLoader/VTK)

		include(${VTK_USE_FILE})
		set (EXTRA_LIBS ${EXTRA_LIBS} ${VTK_LIBRARIES})
	else()
		message(FATAL_ERROR "Could not find vtk!!!") 
	endif()
endif()


# NYX Loader
option(CBENCH_LOADER_NYX "Enables the reading of NyX HDF5 data" OFF)
if (CBENCH_LOADER_NYX)
	set(HDF5_USE_STATIC_LIBRARIES ON)
	find_package(HDF5 REQUIRED COMPONENTS C)

	if(HDF5_FOUND)
		add_definitions (-DCBENCH_HAS_NYX)
		include_directories(${CMAKE_SOURCE_DIR}/dataLoader/NYX)

		if (${CBENCH_PLATFORM} MATCHES "Cori")
			include_directories("/usr/common/software/hdf5-parallel/1.10.1/gnu/include")
			link_directories("/usr/common/software/hdf5-parallel/1.10.1//gnu/lib")
			set (EXTRA_LIBS ${EXTRA_LIBS} dl)
			set (EXTRA_LIBS ${EXTRA_LIBS} m)
			set (EXTRA_LIBS ${EXTRA_LIBS} z)
			set (EXTRA_LIBS ${EXTRA_LIBS} hdf5.a)
			set (EXTRA_LIBS ${EXTRA_LIBS} hdf5_hl.a)
		else()
			include_directories(${HDF5_INCLUDE_DIRS})
			set(EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARIES})  
		endif()
	else()
		message(FATAL_ERROR "HDF5 library not found")
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
	if (${CBENCH_PLATFORM} MATCHES "Cori")
		set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
		set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
	else()
		set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
		set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
		set(CMAKE_C_COMPILER mpicc)
		set(CMAKE_CXX_COMPILER mpicxx)
	endif()
else()
	message(FATAL_ERROR "Could not find any MPI installation!!!") 
endif()



# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



# Build executable
add_executable(CBench main.cpp)
target_link_libraries(CBench ${MPI_LIBRARIES} ${EXTRA_LIBS} GenericIO)

#add_executable(TestAMReX test.cpp)
#target_link_libraries(TestAMReX ${MPI_LIBRARIES} ${EXTRA_LIBS} GenericIO)


# Install instructions
install(TARGETS CBench DESTINATION .)


# GenericIO
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/genericio)
add_definitions(-DGENERICIO_NO_COMPRESSION)
#add_definitions(-DGENERICIO_NO_MPI)	# might need turning on sometimnes
add_library(GenericIO STATIC "${CMAKE_SOURCE_DIR}/thirdparty/genericio/GenericIO.cxx")
set_property(TARGET GenericIO PROPERTY POSITION_INDEPENDENT_CODE ON)
install(TARGETS GenericIO DESTINATION .)


# Prevent in source builds since they are evil
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
	message(FATAL_ERROR "In source build is not supported, please restore modified files, 
						delete CMake generated things (CMakeFiles folder and CMakeCache.txt)
						, and build out of source") 
endif()


# OS things
if(WIN32)
	add_definitions(-DNOMINMAX)
	file(GLOB_RECURSE EXTRA_DLLS
		"${CMAKE_SOURCE_DIR}/../ExternalDependencies/*.dll"
	)
	message( "Config Types: ${CMAKE_CONFIGURATION_TYPES} Libraries used: ${EXTRA_DLLS}")
	foreach(DLL ${EXTRA_DLLS})
		add_custom_command(TARGET CBench POST_BUILD
			COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${DLL}" "${CMAKE_BINARY_DIR}/$<CONFIG>/"
		)
	endforeach()
endif()