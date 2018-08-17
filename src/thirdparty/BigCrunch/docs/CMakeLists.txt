# doxygen CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# project name
project(doxygen)

# find package
find_package(Doxygen REQUIRED)

# set input and output files
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# configure input file
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

# add custom command to start doxygen build
add_custom_target(${PROJECT_NAME} ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating ${PROJECT_NAME} API documentation with Doxygen"
    VERBATIM
)