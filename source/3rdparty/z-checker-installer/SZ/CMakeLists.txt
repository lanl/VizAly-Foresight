# project name
project(sz)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add include directories
target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/sz/include
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/zlib/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/sz/.libs/libsz.a
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/sz/.libs/libsz.so
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/zlib/.libs/libzlib.a
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/zlib/.libs/libzlib.so
)

# install instructions
file(GLOB
	PROJECT_SHARED_LIB_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/sz/.libs/*.so"
	"${CMAKE_CURRENT_SOURCE_DIR}/sz/.libs/*.so.*"
	"${CMAKE_CURRENT_SOURCE_DIR}/zlib/.libs/*.so"
	"${CMAKE_CURRENT_SOURCE_DIR}/zlib/.libs/*.so.*"
)

install(FILES
	${PROJECT_SHARED_LIB_FILES}
	DESTINATION lib/3rdparty/${PROJECT_NAME}
)
