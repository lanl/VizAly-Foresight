
if (CBENCH_ENABLE_LOSSY_WAVE)

    #Asks for LossyWave library
    find_library(LOSSYWAVE_LIBRARY ${LOSSYWAVE_LIBRARY} CACHE STRING "LOSSYWAVE_LIBRARY library")

    if (NOT LOSSYWAVE_LIBRARY)
    	message(FATAL_ERROR "Could NOT find LossyWave")
    else()
		# Find lz4 dependency
		find_library(LOSSYWAVE_LZ4_LIBRARY ${LOSSYWAVE_LZ4_LIBRARY} CACHE STRING "LOSSYWAVE_LZ4_LIBRARY library")

		#Asks for LossyWave include files and includes it in the project
		set(LOSSYWAVE_INCLUDE_PATH ${LOSSYWAVE_INCLUDE_PATH} CACHE PATH "LOSSYWAVE include directory")
		set(LOSSYWAVE_LZ4_INCLUDE_PATH ${LOSSYWAVE_LZ4_INCLUDE_PATH} CACHE PATH "LOSSYWAVE LZ4 include directory")
		set(LOSSYWAVE_INCLUDE_PATH ${LOSSYWAVE_INCLUDE_PATH} ${LOSSYWAVE_LZ4_INCLUDE_PATH})

	    set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_LOSSY_WAVE")
	   	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${LOSSYWAVE_INCLUDE_PATH})
	    set_property(TARGET CBench APPEND PROPERTY linkLibraries ${LOSSYWAVE_LIBRARY} ${LOSSYWAVE_LZ4_LIBRARY})

	endif()
endif()
