
if (CBENCH_ENABLE_SZ)

    #Asks for SZ library + set that in EXTRA_LIBS
    find_library(SZ_LIBRARY ${SZ_LIBRARY} CACHE STRING "SZ library")
    if (NOT SZ_LIBRARY)
    	message(FATAL_ERROR "Could NOT find SZ")
    else()
    	#Find other libraries
    	find_library(ZLIB_LIBRARY ${ZLIB_LIBRARY} CACHE STRING "ZLIB library")
    	find_library(ZSTD_LIBRARY ${ZSTD_LIBRARY} CACHE STRING "SZ library")

    	#Asks for zfp include files + includes it in the project
    	set (SZ_INCLUDE_PATH ${SZ_INCLUDE_PATH} CACHE PATH "SZ include directory")
    	set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_SZ")
    	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${SZ_INCLUDE_PATH})
    	set_property(TARGET CBench APPEND PROPERTY linkLibraries ${SZ_LIBRARY} ${ZLIB_LIBRARY} ${ZSTD_LIBRARY})
    endif()
endif()