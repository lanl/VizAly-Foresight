
if (CBENCH_ENABLE_SZ)

    #Asks for SZ library + set that in EXTRA_LIBS
    find_library(SZ_LIBRARY ${SZ_LIBRARY} CACHE STRING "SZ library")

	if (SZ_LIBRARY)
		#Find other libraries
		find_library(ZLIB_LIBRARY ${ZLIB_LIBRARY} CACHE STRING "ZLIB library")
		find_library(ZSTD_LIBRARY ${ZSTD_LIBRARY} CACHE STRING "ZSTD library")

		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_SZ)
		target_include_directories(CBench PRIVATE ${SZ_INCLUDE_PATH})
		target_include_directories(GenericIO PRIVATE ${SZ_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY LIB "${SZ_LIBRARY}")
		set_property(TARGET CBench APPEND PROPERTY LIB "${ZLIB_LIBRARY}")
		set_property(TARGET CBench APPEND PROPERTY LIB "${ZSTD_LIBRARY}")
	else()
		message(FATAL_ERROR "Could NOT find SZ")
	endif()
endif()