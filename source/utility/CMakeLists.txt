# library CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.5.1)

# project name
project(utility)

# capital project name string
string(TOUPPER ${PROJECT_NAME} NAME_ID)

# option for STATIC or SHARED
option(${NAME_ID}_BUILD_STATIC
	"Build a static version of this library."
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
	src/benchmark.cpp
	src/error.cpp
	src/fio.cpp
	src/stopwatch.cpp
)

# add header files
set(PROJECT_HEADER_FILES
)

# add interface header files
set(PROJECT_INTERFACE_HEADER_FILES
	include/benchmark.hpp
	include/error.hpp
	include/fio.hpp
	include/stopwatch.hpp
	include/log.hpp
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
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	PUBLIC mpi
)

# generate visibility header
#generate_export_header(${PROJECT_NAME})

# install instructions
if(${NAME_ID}_BUILD_STATIC AND ${NAME_ID}_INSTALL)
	install(TARGETS ${PROJECT_NAME}
		DESTINATION lib
	)
else()
	if(WIN32)
		install(TARGETS ${PROJECT_NAME}
			DESTINATION bin
		)
	else()
		install(TARGETS ${PROJECT_NAME}
			DESTINATION lib
		)
	endif()
endif()

if(${NAME_ID}_INSTALL)
	install(FILES
		${PROJECT_INTERFACE_HEADER_FILES}
		DESTINATION include/${PROJECT_NAME}
	)
endif()
