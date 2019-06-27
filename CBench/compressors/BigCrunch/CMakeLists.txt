
if (CBENCH_ENABLE_BIG_CRUNCH)

    # Asks for BigCrunch library
    find_library(BIGCRUNCH_LIBRARY ${BIGCRUNCH_LIBRARY} CACHE STRING "BIGCRUNCH_LIBRARY library")

    if (BIGCRUNCH_LIBRARY)
	    #Asks for BigCrunch include files and includes it in the project
	    #set (BIGCRUNCH_INCLUDE_PATH ${BIGCRUNCH_INCLUDE_PATH} CACHE PATH "BIGCRUNCH include directory")
		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_BIG_CRUNCH)
		target_include_directories(CBench PRIVATE ${BIGCRUNCH_INCLUDE_PATH})
		target_include_directories(GenericIO PRIVATE ${BIGCRUNCH_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY LIB "${BIGCRUNCH_LIBRARY}")
	else()
		message(FATAL_ERROR "Could NOT find BigCrunch")
	endif()
endif()
