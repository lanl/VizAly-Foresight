option(CBENCH_COMPRESSOR_ZFP "use the zfp compressor" OFF)

if(CBENCH_COMPRESSOR_ZFP)
	set(MY_LIBS "")

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(ZFP_LIBRARY ${ZFP_LIBRARY} CACHE STRING "ZFP_LIBRARY library")

    if (NOT ZFP_LIBRARY)
    	message(FATAL_ERROR "Could NOT find ZFP")
    else()
    	set(MY_LIBS ${MY_LIBS} ${ZFP_LIBRARY})

    	#Asks for zfp include files + includes it in the project
    	set (ZFP_INCLUDE_PATH ${ZFP_INCLUDE_PATH} CACHE PATH "ZFP include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_ZFP")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${ZFP_INCLUDE_PATH})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()