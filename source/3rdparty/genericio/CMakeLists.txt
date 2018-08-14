# library CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.5.1)

# project name
project(genericio)

# captial project name string
string(TOUPPER ${PROJECT_NAME} NAME_ID)

# option for STATIC or SHARED
option(${NAME_ID}_BUILD_STATIC
	"Build a static version of this library"
	ON
)

# option for installing the target
option(${NAME_ID}_INSTALL
	"Install the include headers and the static target (if build) of ${PROJECT_NAME}."
	ON
)

# include modules
#include(GenerateExportHeader)

# add source files
set(PROJECT_SOURCE_FILES
	DataFilter.cxx
	GenericIO.cxx
	GenericIOBase.cxx
	GenericIOMPIReader.cxx
	GenericIOPosixReader.cxx
	GenericIOReader.cxx
	GenericIOUtilities.cxx
)

# add header files
set(PROJECT_HEADER_FILES
#	include/source.h
)

# add interface header files
set(PROJECT_INTERFACE_HEADER_FILES
	include/AssocVector.h
	include/CRC64.h
	include/DataFilter.h
	include/GenericIOBase.h
	include/GenericIOCompress.h
	include/GenericIODefinitions.hpp
	include/GenericIO.h
	include/GenericIOMPIReader.h
	include/GenericIOPosixReader.h
	include/GenericIOReader.h
	include/GenericIOUtilities.h
)

# add library targets
if(${NAME_ID}_BUILD_STATIC)
	add_library(${PROJECT_NAME} STATIC
		${PROJECT_SOURCE_FILES}
		${PROJECT_HEADER_FILES}
		${PROJECT_INTERFACE_HEADER_FILES}
	)
else()
	add_library(${PROJECT_NAME} SHARED
		${PROJECT_SOURCE_FILES}
		${PROJECT_HEADER_FILES}
		${PROJECT_INTERFACE_HEADER_FILES}
	)
endif()

# set c++ standard
set_target_properties(${PROJECT_NAME} PROPERTIES
	CXX_STANDARD 11
	CXX_STANDARD_REQUIRED ON
	CXX_EXTENSIONS OFF
)

# add include directories
target_include_directories(${PROJECT_NAME}
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	PUBLIC mpi
	PUBLIC openmp
)

# add compile definitions
# if no definitions need to be set, this command needs to be commented out
target_compile_definitions(${PROJECT_NAME}
	PUBLIC __STDC_CONSTANT_MACROS
)

# add compile flags
# if no flags need to be set, this command needs to be commented out
#target_compile_options(${PROJECT_NAME}
#)

# generate visibility header
#generate_export_header(${PROJECT_NAME})

# install instructions
if(${NAME_ID}_BUILD_STATIC AND ${NAME_ID}_INSTALL)
	install(TARGETS ${PROJECT_NAME}
		DESTINATION lib/3rdparty/${PROJECT_NAME}
	)
elseif(NOT ${NAME_ID}_BUILD_STATIC)
	if(WIN32)
		install(TARGETS ${PROJECT_NAME}
			DESTINATION bin
		)
	else()
		install(TARGETS ${PROJECT_NAME}
			DESTINATION lib/3rdparty/${PROJECT_NAME}
		)
	endif()
endif()

if(${NAME_ID}_INSTALL)
	install(FILES
		${PROJECT_INTERFACE_HEADER_FILES}
		DESTINATION include/3rdparty/${PROJECT_NAME}
	)
endif()
