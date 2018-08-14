# project name
project(mpi)

# use old policy because of bug in cmake
cmake_policy(SET CMP0004 OLD)

# find package
find_package(MPI REQUIRED)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add compile definitions
# if no definitions need to be set, this command needs to be commented out
#target_compile_definitions(${PROJECT_NAME}
#)

# add compile flags
# if no flags need to be set, this command needs to be commented out
target_compile_options(${PROJECT_NAME}
	INTERFACE ${MPI_CXX_COMPILE_FLAGS}
)

# add include directories
target_include_directories(${PROJECT_NAME}
	INTERFACE ${MPI_CXX_INCLUDE_PATH}
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	INTERFACE ${MPI_CXX_LINK_FLAGS}
	INTERFACE ${MPI_CXX_LIBRARIES}
)
