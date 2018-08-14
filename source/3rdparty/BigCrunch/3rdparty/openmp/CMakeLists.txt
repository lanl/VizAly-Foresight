#external package CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# project name
project(openmp LANGUAGES CXX)

# find package
find_package(OpenMP REQUIRED)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add compile flags
target_compile_options(${PROJECT_NAME}
    INTERFACE ${OpenMP_CXX_FLAGS}
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
    INTERFACE ${OpenMP_CXX_LIBRARIES}
)