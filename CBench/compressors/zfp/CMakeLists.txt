
if(CBENCH_ENABLE_ZFP)

    #Asks for zfp library + set that in EXTRA_LIBS
    find_library(ZFP_LIBRARY ${ZFP_LIBRARY} CACHE STRING "ZFP_LIBRARY library")

    if (NOT ZFP_LIBRARY)
    	message(FATAL_ERROR "Could NOT find ZFP")
    else()
    	#Asks for zfp include files + includes it in the project
    	set(ZFP_INCLUDE_PATH ${ZFP_INCLUDE_PATH} CACHE PATH "ZFP include directory")
    	set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_ZFP")
    	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${ZFP_INCLUDE_PATH})
    	set_property(TARGET CBench APPEND PROPERTY linkLibraries ${MY_LIBS} ${ZFP_LIBRARY})
    endif()
endif()