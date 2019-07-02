
if(CBENCH_ENABLE_ZFP)

    find_library(ZFP_LIBRARY ${ZFP_LIBRARY} CACHE STRING "ZFP_LIBRARY library")

	if (ZFP_LIBRARY)
		set(ZFP_INCLUDE_PATH ${ZFP_INCLUDE_PATH} CACHE PATH "ZFP include directory")
		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_ZFP)
		target_include_directories(CBench PRIVATE ${ZFP_INCLUDE_PATH})
		target_include_directories(GenericIO PRIVATE ${ZFP_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY LIB "${ZFP_LIBRARY}")
	else()
		message(FATAL_ERROR "Could NOT find ZFP")
	endif()
endif()