
if (CBENCH_ENABLE_MGARD)

    #Asks for MGARD library
    find_library(MGARD_LIBRARY ${MGARD_LIBRARY} CACHE STRING "MGARD_LIBRARY library")

    if (NOT MGARD_LIBRARY)
    	message(FATAL_ERROR "Could NOT find MGARD")
    else()
	    #Asks for MGARD include files and includes it in the project
	    set (MGARD_INCLUDE_PATH ${MGARD_INCLUDE_PATH} CACHE PATH "MGARD include directory")
	    set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_MGARD")
	   	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${MGARD_INCLUDE_PATH})
	    set_property(TARGET CBench APPEND PROPERTY linkLibraries ${MGARD_LIBRARY})
	endif()
endif()
