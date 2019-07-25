option(CBENCH_ENABLE_SNWPAC "Enable snwpac." OFF)

if (CBENCH_ENABLE_SNWPAC)
	set(MY_LIBS "")

    #Asks for SNWPAC library
    find_library(SNWPAC_LIBRARY ${SNWPAC_LIBRARY} CACHE STRING "SNWPAC_LIBRARY library")

    if (NOT SNWPAC_LIBRARY)
    	message(FATAL_ERROR "Could NOT find SNWPAC")
    else()
            
    	set(MY_LIBS ${MY_LIBS} ${SNWPAC_LIBRARY})
	    set(MY_LIBS ${MY_LIBS} -fopenmp)
       
	    #Asks for SNWPAC include files and includes it in the project
	    set (SNWPAC_INCLUDE_PATH ${SNWPAC_INCLUDE_PATH} CACHE PATH "SNWPAC include directory")

		# Find lz4 dependency (only on Win32)
		#if(WIN32)
			find_library(SNWPAC_LZ4_LIBRARY ${SNWPAC_LZ4_LIBRARY} CACHE STRING "SNWPAC_LZ4_LIBRARY library")
			set(MY_LIBS ${MY_LIBS} ${SNWPAC_LZ4_LIBRARY})
			set(SNWPAC_LZ4_INCLUDE_PATH ${SNWPAC_LZ4_INCLUDE_PATH} CACHE PATH "SNWPAC LZ4 include directory")
			set(SNWPAC_INCLUDE_PATH ${SNWPAC_INCLUDE_PATH} ${SNWPAC_LZ4_INCLUDE_PATH})
		#endif()

	    set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_SNWPAC")
	   	set_property(GLOBAL APPEND PROPERTY includeDirectories ${SNWPAC_INCLUDE_PATH})
	    set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
        
	endif()
endif()
