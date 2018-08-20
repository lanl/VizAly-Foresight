# library CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# project name
project(bigcrunch LANGUAGES CXX)

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

# option for building tests
# option(${NAME_ID}_TESTS
# 	"Build unit tests for ${PROJECT_NAME}."
# 	OFF
# )

# option for building doxygen documentation
option(${NAME_ID}_DOXYGEN
	"Build doxygen documentation for ${PROJECT_NAME}."
	OFF
)

# option for building doxygen documentation only in release mode
option(${NAME_ID}_DOXYGEN_RELEASE_ONLY
	"Build doxygen documentation for ${PROJECT_NAME} only in Release mode."
	OFF
)

# add source files
set(PROJECT_SOURCE_FILES
	src/bigcrunch.cpp
	src/bit_ops.cpp
	src/blosc_wrapper.cpp
	src/byte_buffer.cpp
	src/darray.cpp
	src/fpr_data.cpp
	src/fpr_stats.cpp
	src/fpr.cpp
)

# add headder files
set(PROJECT_HEADER_FILES
	include/bigcrunch/bit_ops.hpp
	include/bigcrunch/blosc_wrapper.hpp
	include/bigcrunch/byte_buffer.hpp
	include/bigcrunch/fpr_data.hpp
	include/bigcrunch/fpr_decoder.hpp
	include/bigcrunch/fpr_encoder.hpp
	include/bigcrunch/fpr_ops.hpp
	include/bigcrunch/fpr_stats.hpp
	include/bigcrunch/fpr.hpp
	include/bigcrunch/detail/byte_buffer.inl
	include/bigcrunch/detail/fpr_decoder.inl
	include/bigcrunch/detail/fpr_encoder.inl
	include/bigcrunch/detail/fpr_ops.inl
)

# add interface header files
set(PROJECT_INTERFACE_HEADER_FILES
	include/bigcrunch/bigcrunch.hpp
	include/bigcrunch/darray.hpp
	include/bigcrunch/types.hpp
	include/bigcrunch/detail/darray.inl
)

# 3rdparty libraries
add_subdirectory(3rdparty/openmp)
if (NOT TARGET blosc)
    add_subdirectory(3rdparty/blosc)
endif() 

# test targets
# if(${NAME_ID}_TESTS)
add_subdirectory(tests)
# endif()

# add doxygen target
if(${NAME_ID}_DOXYGEN)
	if(${NAME_ID}_DOXYGEN_RELEASE_ONLY)
		if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
			add_subdirectory(docs)
		endif()
	else()
		add_subdirectory(docs)
	endif()
endif()

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
	PUBLIC openmp
	PUBLIC blosc
)

# install instructions
if(${NAME_ID}_BUILD_STATIC AND ${NAME_ID}_INSTALL)
	install(TARGETS ${PROJECT_NAME}
		DESTINATION lib
	)
elseif(NOT ${NAME_ID}_BUILD_STATIC)
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
