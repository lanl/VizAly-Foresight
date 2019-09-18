option(CBENCH_COMPRESSOR_SZ "Enable SZ." OFF)

if (CBENCH_COMPRESSOR_SZ)
	set(MY_LIBS "")

    #Asks for SZ library + set that in EXTRA_LIBS
    find_library(SZ_LIBRARY ${SZ_LIBRARY} CACHE STRING "SZ library")
    if (NOT SZ_LIBRARY)
    	message(FATAL_ERROR "Could NOT find SZ")
    else()
    	set(MY_LIBS ${MY_LIBS} ${SZ_LIBRARY})

    	#Find other libraries
    	find_library(ZLIB_LIBRARY ${ZLIB_LIBRARY} CACHE STRING "ZLIB library")
    	set(MY_LIBS ${MY_LIBS} ${ZLIB_LIBRARY})

    	find_library(ZSTD_LIBRARY ${ZSTD_LIBRARY} CACHE STRING "SZ library")
    	set(MY_LIBS ${MY_LIBS} ${ZSTD_LIBRARY})

    	#Asks for zfp include files + includes it in the project
    	set (SZ_INCLUDE_PATH ${SZ_INCLUDE_PATH} CACHE PATH "SZ include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_SZ")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${SZ_INCLUDE_PATH})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()