# external project CMakeLists.txt template

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# project name
project(blosc)

# enable external project cmake module
include(ExternalProject)

# download external project
ExternalProject_Add(
    ${PROJECT_NAME}-external
    #URL https://github.com/Blosc/c-blosc/archive/master.zip
	GIT_REPOSITORY https://github.com/Blosc/c-blosc.git
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-external
    # disable install step
    INSTALL_COMMAND ""
)

# get external project source and binary directories
ExternalProject_Get_Property(${PROJECT_NAME}-external
    SOURCE_DIR
    BINARY_DIR
)

# create source directory with include path.
# this is necessary to configure interface include directories.
execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SOURCE_DIR}/blosc
)

# create imported target to be used as a dependency by programs
add_library(${PROJECT_NAME} STATIC IMPORTED GLOBAL)

# set project properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    IMPORTED_LOCATION ${BINARY_DIR}/blosc/libblosc.a
)

# add external project as dependency
add_dependencies(${PROJECT_NAME}
    ${PROJECT_NAME}-external
)
# target_link_libraries(${PROJECT_NAME}
#     INTERFACE ${PROJECT_NAME}-external
# )

# add external project include directories
target_include_directories(${PROJECT_NAME}
    INTERFACE ${SOURCE_DIR}/blosc
)