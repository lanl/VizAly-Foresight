
if (CBENCH_ENABLE_ISABELA)

    #Asks for Isabela library + set that in EXTRA_LIBS
    find_library(ISABELA_LIBRARY ${ISABELA_LIBRARY} CACHE STRING "ISABELA library")
    if (NOT ISABELA_LIBRARY)
    	message(FATAL_ERROR "Could NOT find ISABELA")
    else()
		#Find other libraries
		find_package(GSL REQUIRED)

    	#Asks for ISABELA include files + includes it in the project
    	set (ISABELA_INCLUDE_PATH ${ISABELA_INCLUDE_PATH} CACHE PATH "ISABELA include directory")

    	set_property(TARGET CBench APPEND PROPERTY definitions "-DCBENCH_HAS_ISABELA")
    	set_property(TARGET CBench APPEND PROPERTY includeDirectories ${ISABELA_INCLUDE_PATH})
		set_property(TARGET CBench APPEND PROPERTY includeDirectories ${GSL_INCLUDE_DIRS})
    	set_property(TARGET CBench APPEND PROPERTY linkLibraries ${ISABELA_LIBRARY} ${GSL_LIBRARIES})
    endif()
endif()