
if (CBENCH_ENABLE_FPZIP)
    #Asks for fpzip library + set that in EXTRA_LIBS
    find_library(FPZIP_LIBRARY ${FPZIP_LIBRARY} CACHE STRING "FPZIP library")
    if (NOT FPZIP_LIBRARY)
    	message(FATAL_ERROR "Could NOT find FPZIP")
    else()
    	#Asks for fpzip include files + includes it in the project
    	set (FPZIP_INCLUDE_PATH ${FPZIP_INCLUDE_PATH} CACHE PATH "FPZIP include directory")

    	set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_FPZIP")
    	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${FPZIP_INCLUDE_PATH})
    	set_property(TARGET CBench APPEND PROPERTY linkLibraries ${FPZIP_LIBRARY})
    endif()
endif()