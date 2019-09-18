option(CBENCH_COMPRESSOR_BLOSC "Enable BLOSC." OFF)

if(CBENCH_COMPRESSOR_BLOSC)
	set(MY_LIBS "")

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(BLOSC_LIBRARY ${BLOSC_LIBRARY} CACHE STRING "BLOSC_LIBRARY library")

    if (NOT BLOSC_LIBRARY)
    	message(FATAL_ERROR "Could NOT find BLOSC")
    else()
    	set(MY_LIBS ${MY_LIBS} ${BLOSC_LIBRARY})

    	#Asks for zfp include files + includes it in the project
    	set (BLOSC_INCLUDE_PATH ${BLOSC_INCLUDE_PATH} CACHE PATH "BLOSC include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_BLOSC")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${BLOSC_INCLUDE_PATH})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()