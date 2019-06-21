
if (CBENCH_ENABLE_MGARD)

    find_library(MGARD_LIBRARY ${MGARD_LIBRARY} CACHE STRING "MGARD_LIBRARY library")

	if (MGARD_LIBRARY)
		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_MGARD)
		target_include_directories(CBench PRIVATE ${MGARD_INCLUDE_PATH})
		target_include_directories(GenericIO PRIVATE ${MGARD_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY LIB "${MGARD_LIBRARY}")
	else()
		message(FATAL_ERROR "Could NOT find MGARD")
	endif()
endif()
