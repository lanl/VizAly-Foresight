option(CBENCH_ENABLE_MGARD "Enable MGARD." OFF)

if (CBENCH_ENABLE_MGARD)
	set(MY_LIBS "")

    #Asks for MGARD library
    find_library(MGARD_LIBRARY ${MGARD_LIBRARY} CACHE STRING "MGARD_LIBRARY library")

    if (NOT MGARD_LIBRARY)
    	message(FATAL_ERROR "Could NOT find MGARD")
    else()
            
    	set(MY_LIBS ${MY_LIBS} ${MGARD_LIBRARY})
	    set(MY_LIBS ${MY_LIBS} -fopenmp)
       
	    #Asks for MGARD include files and includes it in the project
	    set (MGARD_INCLUDE_PATH ${MGARD_INCLUDE_PATH} CACHE PATH "MGARD include directory")

	    set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_MGARD")
	   	set_property(GLOBAL APPEND PROPERTY includeDirectories ${MGARD_INCLUDE_PATH})
	    set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
        
	endif()
endif()
