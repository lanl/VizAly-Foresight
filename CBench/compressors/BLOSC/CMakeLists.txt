
if(CBENCH_ENABLE_BLOSC)
    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(BLOSC_LIBRARY ${BLOSC_LIBRARY} CACHE STRING "BLOSC_LIBRARY library")

    if (NOT BLOSC_LIBRARY)
    	message(FATAL_ERROR "Could NOT find BLOSC")
    else()
    	#Asks for zfp include files + includes it in the project
    	set (BLOSC_INCLUDE_PATH ${BLOSC_INCLUDE_PATH} CACHE PATH "BLOSC include directory")

    	set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_BLOSC")
    	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${BLOSC_INCLUDE_PATH})
    	set_property(TARGET CBench APPEND PROPERTY linkLibraries ${BLOSC_LIBRARY})
    endif()
endif()