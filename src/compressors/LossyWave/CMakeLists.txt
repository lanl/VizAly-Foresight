option(CBENCH_ENABLE_LOSSY_WAVE "Enable LossyWave." ON)

if (CBENCH_ENABLE_LOSSY_WAVE)
	set(MY_LIBS "")

    #Asks for BigCrucnch library
    find_library(LOSSYWAVE_LIBRARY ${LOSSYWAVE_LIBRARY} CACHE STRING "LOSSYWAVE_LIBRARY library")

    if (NOT LOSSYWAVE_LIBRARY)
    	message(FATAL_ERROR "Could NOT find LossyWave")
    else()
    	set(MY_LIBS ${MY_LIBS} ${LOSSYWAVE_LIBRARY})
	    set(MY_LIBS ${MY_LIBS} -fopenmp)
	   
	    #Asks for LossyWave include files and includes it in the project
	    set (LOSSYWAVE_INCLUDE_PATH ${LOSSYWAVE_INCLUDE_PATH} CACHE STRING "LOSSYWAVE include directory")


	    set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_LOSSY_WAVE")
	   	set_property(GLOBAL APPEND PROPERTY includeDirectories ${LOSSYWAVE_INCLUDE_PATH})
	    set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
	endif()
endif()
