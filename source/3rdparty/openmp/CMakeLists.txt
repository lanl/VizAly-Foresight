# project name
project(openmp)

# find packages
find_package(OpenMP REQUIRED)

message(${OpenMP_CXX_FLAGS})

# add header-only library
add_library(${PROJECT_NAME} INTERFACE)

# add compile flags
target_compile_options(${PROJECT_NAME}
	INTERFACE ${OpenMP_CXX_FLAGS}
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	INTERFACE ${OpenMP_CXX_FLAGS}
)
