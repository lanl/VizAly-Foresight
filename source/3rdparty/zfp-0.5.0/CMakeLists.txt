# project name
project(zfp)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add include directories
target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

# add libary dependencies
target_link_libraries(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/lib/libzfp.a
)
