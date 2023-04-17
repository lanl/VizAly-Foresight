option(CBENCH_COMPRESSOR_SZ3 "Enable SZ3." OFF)

if (CBENCH_COMPRESSOR_SZ3)
	set(MY_LIBS "")

    #Asks for SZ library + set that in EXTRA_LIBS
    find_library(SZ3_LIBRARY ${SZ3_LIBRARY} CACHE STRING "SZ3 library")
    if (NOT SZ3_LIBRARY)
    	message(FATAL_ERROR "Could NOT find SZ3")
    else()
    	set(MY_LIBS ${MY_LIBS} ${SZ3_LIBRARY})

    	#Find other libraries
    	#find_library(ZLIB_LIBRARY ${ZLIB_LIBRARY} CACHE STRING "ZLIB library")
    	#set(MY_LIBS ${MY_LIBS} ${ZLIB_LIBRARY})

		#find_library(ZSTD_LIBRARY ${ZSTD_LIBRARY} CACHE STRING "SZ library")
    	#set(MY_LIBS ${MY_LIBS} ${ZSTD_LIBRARY})

    	#Asks for zfp include files + includes it in the project
    	set (SZ3_INCLUDE_PATH ${SZ3_INCLUDE_PATH} CACHE PATH "SZ3 include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_SZ3")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${SZ3_INCLUDE_PATH})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()
