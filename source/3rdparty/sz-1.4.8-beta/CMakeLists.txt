# project name
project(sz)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add include directories
target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/build/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/libsz.a
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/libsz.so
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/libzlib.a
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/libzlib.so
)

# install instructions
file(GLOB
	PROJECT_SHARED_LIB_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/build/lib/*.so"
	"${CMAKE_CURRENT_SOURCE_DIR}/build/lib/*.so.*"
)

install(FILES
	${PROJECT_SHARED_LIB_FILES}
	DESTINATION lib/3rdparty/${PROJECT_NAME}
)
