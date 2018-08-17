#external package CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.10)

# project name
project(thrust LANGUAGES CXX CUDA)

# option for thrust device backends
set(THRUST_DEVICE_SYSTEM "OpenMP" CACHE STRING "Option for setting the backend used by thrust.")
set_property(CACHE THRUST_DEVICE_SYSTEM PROPERTY STRINGS OpenMP Cuda)

# find package
find_package(OpenMP REQUIRED)

# add header-only target
add_library(${PROJECT_NAME} INTERFACE)

# add compile flags
if(${THRUST_DEVICE_SYSTEM} STREQUAL "OpenMP")
    target_compile_options(${PROJECT_NAME}
        INTERFACE -Xcompiler "${OpenMP_CXX_FLAGS}"
        INTERFACE -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
    )
elseif(${THRUST_DEVICE_SYSTEM} STREQUAL "Cuda")
    target_compile_options(${PROJECT_NAME}
        INTERFACE -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA
    )
else()
    message(FATAL_ERROR "Unsupported thrust device backend!")
endif()

target_compile_options(${PROJECT_NAME}
    INTERFACE -DCUDA_SEPARABLE_COMPILATION=ON
)

if(APPLE)
    target_compile_options(${PROJECT_NAME}
        INTERFACE -DBUILD_RPATH=${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}
    )
endif()

# add library dependencies
target_link_libraries(${PROJECT_NAME}
    INTERFACE ${OpenMP_CXX_LIBRARIES}
)