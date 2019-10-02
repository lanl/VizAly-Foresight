option(CBENCH_COMPRESSOR_SZ_GPU "Enable SZ_gpu." OFF)

if (CBENCH_COMPRESSOR_SZ_GPU)
	set(MY_LIBS "")

    #Asks for SZ_gpu library + set that in EXTRA_LIBS
    find_library(SZ_GPU_LIBRARY ${SZ_GPU_LIBRARY} CACHE STRING "SZ_GPU library")
    if (NOT SZ_GPU_LIBRARY)
    	message(FATAL_ERROR "Could NOT find SZ_GPU")
    else()
    	set(MY_LIBS ${MY_LIBS} ${SZ_GPU_LIBRARY})

    	#Find other libraries
    	find_library(ZLIB_GPU_LIBRARY ${ZLIB_GPU_LIBRARY} CACHE STRING "ZLIB_GPU library")
    	set(MY_LIBS ${MY_LIBS} ${ZLIB_GPU_LIBRARY})

    	find_library(ZSTD_GPU_LIBRARY ${ZSTD_GPU_LIBRARY} CACHE STRING "SZ_GPU library")
    	set(MY_LIBS ${MY_LIBS} ${ZSTD_GPU_LIBRARY})

    	#Asks for zfp include files + includes it in the project
    	set (SZGPU_INCLUDE_PATH ${SZ_GPU_INCLUDE_PATH} CACHE PATH "SZ_GPU include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_SZ_GPU")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${SZ_GPU_INCLUDE_PATH})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()