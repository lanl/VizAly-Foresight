# project name
project(fpzip)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add include directories
target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/lib/libfpzip.a
)
