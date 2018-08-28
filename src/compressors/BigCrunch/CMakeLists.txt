option(CBENCH_ENABLE_BIG_CRUNCH "Enable BigCrunch." OFF)

if (CBENCH_ENABLE_BIG_CRUNCH)
	set(MY_LIBS "")

    #Asks for BigCrucnch library
    find_library(BIGCRUNCH_LIBRARY ${BIGCRUNCH_LIBRARY} CACHE STRING "BIGCRUNCH_LIBRARY library")

    if (NOT BIGCRUNCH_LIBRARY)
    	message(FATAL_ERROR "Could NOT find BigCrunch")
    else()
    	set(MY_LIBS ${MY_LIBS} ${BIGCRUNCH_LIBRARY})

	    find_library(BIGCRUNCH_BLOSC_LIBRARY ${BIGCRUNCH_BLOSC_LIBRARY} CACHE STRING "BIGCRUNCH_LIBRARY library")
	    set(MY_LIBS ${MY_LIBS} ${BIGCRUNCH_BLOSC_LIBRARY})

	    set(MY_LIBS ${MY_LIBS} -fopenmp)

	   
	    #Asks for BigCrunch include files and includes it in the project
	    set (BIGCRUNCH_INCLUDE_PATH ${BIGCRUNCH_INCLUDE_PATH} CACHE STRING "BIGCRUNCH include directory")


	    set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_BIG_CRUNCH")
	   	set_property(GLOBAL APPEND PROPERTY includeDirectories ${BIGCRUNCH_INCLUDE_PATH})
	    set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
	endif()
endif()