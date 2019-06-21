
if (CBENCH_ENABLE_LOSSY_WAVE)

    #Asks for LossyWave library
    find_library(LOSSYWAVE_LIBRARY ${LOSSYWAVE_LIBRARY} CACHE STRING "LOSSYWAVE_LIBRARY library")

	if (LOSSYWAVE_LIBRARY)
		# Find lz4 dependency
		find_library(LOSSYWAVE_LZ4_LIBRARY ${LOSSYWAVE_LZ4_LIBRARY} CACHE STRING "LOSSYWAVE_LZ4_LIBRARY library")

		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_LOSSY_WAVE)
		target_include_directories(CBench PRIVATE ${LOSSYWAVE_INCLUDE_PATH} ${LOSSYWAVE_LZ4_INCLUDE_PATH})
		target_include_directories(GenericIO PRIVATE ${LOSSYWAVE_INCLUDE_PATH} ${LOSSYWAVE_LZ4_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY LIB "${LOSSYWAVE_LIBRARY}")
		set_property(TARGET CBench APPEND PROPERTY LIB "${LOSSYWAVE_LZ4_LIBRARY}")
	else()
		message(FATAL_ERROR "Could NOT find LossyWave")
	endif()
endif()
