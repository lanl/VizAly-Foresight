
if (CBENCH_ENABLE_ISABELA)

    find_library(ISABELA_LIBRARY ${ISABELA_LIBRARY} CACHE STRING "ISABELA library")

	if (ISABELA_LIBRARY)
		#Find other libraries
		find_package(GSL REQUIRED)
		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_ISABELA)
		target_include_directories(CBench PRIVATE ${ISABELA_INCLUDE_PATH} ${GSL_INCLUDE_DIRS})
		target_include_directories(GenericIO PRIVATE ${ISABELA_INCLUDE_PATH} ${GSL_INCLUDE_DIRS})
		set_property(TARGET CBench APPEND PROPERTY LIB "${ISABELA_LIBRARY}")
		set_property(TARGET CBench APPEND PROPERTY LIB "${GSL_LIBRARIES}")
	else()
		message(FATAL_ERROR "Could NOT find ISABELA")
	endif()
endif()