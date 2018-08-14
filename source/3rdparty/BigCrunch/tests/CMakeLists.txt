# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6)

# option for building unit tests
option(${NAME_ID}_UNIT_TESTS
    "Build unit tests for ${PROJECT_NAME}."
    OFF
)

# option for building integration test
option(${NAME_ID}_INTEGRATION_TESTS
    "Build integration tests for ${PROJECT_NAME}."
    OFF
)

# add unit test target
if(${NAME_ID}_UNIT_TESTS)
    add_subdirectory(unit)
endif()

# add integration test target
if(${NAME_ID}_INTEGRATION_TESTS)
    add_subdirectory(integration)
endif()