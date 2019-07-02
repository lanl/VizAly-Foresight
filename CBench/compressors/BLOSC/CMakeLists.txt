
if(CBENCH_ENABLE_BLOSC)
    find_library(BLOSC_LIBRARY ${BLOSC_LIBRARY} CACHE STRING "BLOSC_LIBRARY library")

	if (BLOSC_LIBRARY)
		set(BLOSC_INCLUDE_PATH ${BLOSC_INCLUDE_PATH} CACHE PATH "BLOSC include directory")
		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_BLOSC)
		target_include_directories(CBench PRIVATE ${BLOSC_INCLUDE_PATH})
		target_include_directories(GenericIO PRIVATE ${BLOSC_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY LIB "${BLOSC_LIBRARY}")
	else()
		message(FATAL_ERROR "Could NOT find BLOSC")
    endif()
endif()