
if (CBENCH_ENABLE_BIG_CRUNCH)

    # Asks for BigCrunch library
    find_library(BIGCRUNCH_LIBRARY ${BIGCRUNCH_LIBRARY} CACHE STRING "BIGCRUNCH_LIBRARY library")

    if (NOT BIGCRUNCH_LIBRARY)
    	message(FATAL_ERROR "Could NOT find BigCrunch")
    else()
	    #Asks for BigCrunch include files and includes it in the project
	    set (BIGCRUNCH_INCLUDE_PATH ${BIGCRUNCH_INCLUDE_PATH} CACHE PATH "BIGCRUNCH include directory")

		set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_BIG_CRUNCH")
	   	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${BIGCRUNCH_INCLUDE_PATH})
	    set_property(TARGET CBench APPEND PROPERTY linkLibraries ${BIGCRUNCH_LIBRARY})
	endif()
endif()
