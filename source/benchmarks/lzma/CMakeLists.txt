# executable CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.5.1)

# project name
project(lzma_benchmark)

# add source files
set(PROJECT_SOURCE_FILES
	src/compression_benchmark.cpp
	src/decompression_benchmark.cpp
	src/main.cpp
)

# add header files
set(PROJECT_HEADER_FILES
	include/compression_benchmark.hpp
	include/decompression_benchmark.hpp
)

# add executable target
add_executable(${PROJECT_NAME}
	${PROJECT_SOURCE_FILES}
	${PROJECT_HEADER_FILES}
)

# set c++ standard
set_target_properties(${PROJECT_NAME} PROPERTIES
	CXX_STANDARD 11
	CXX_STANDARD_REQUIRED ON
	CXX_EXTENSIONS OFF
)

# add include directories
target_include_directories(${PROJECT_NAME}
	PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	PRIVATE data_loader
	PRIVATE lzma
	PRIVATE utility
)

# install instructions
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)

#install(FILES
#	${OTHER_FILES}
#	DESTINATION somepath
#)
