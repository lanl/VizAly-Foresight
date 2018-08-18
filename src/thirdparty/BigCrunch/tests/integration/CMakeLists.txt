# executable CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# project name
project(integration_test CXX)

# add source files
set(PROJECT_SOURCE_FILES
    src/main.cpp
)

# add header files
set(PROJECT_HEADER_FILES
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
    PRIVATE bigcrunch
)

# install instructions
install(TARGETS ${PROJECT_NAME}
    DESTINATION tests
)