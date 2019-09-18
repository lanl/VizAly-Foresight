option(CBENCH_COMPRESSOR_ZFP_GPU "use the zfp_gpu compressor" OFF)

if(CBENCH_COMPRESSOR_ZFP_GPU)
	set(MY_LIBS "")

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(ZFP_GPU_LIBRARY ${ZFP_GPU_LIBRARY} CACHE STRING "ZFP_GPU_LIBRARY library")

    if (NOT ZFP_GPU_LIBRARY)
    	message(FATAL_ERROR "Could NOT find ZFP_GPU")
    else()
    	set(MY_LIBS ${MY_LIBS} ${ZFP_GPU_LIBRARY})

    	#Asks for zfp include files + includes it in the project
    	set (ZFP_GPU_INCLUDE_PATH ${ZFP_GPU_INCLUDE_PATH} CACHE PATH "ZFP_gpu include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_ZFP_GPU")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${ZFP_GPU_INCLUDE_PATH})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()