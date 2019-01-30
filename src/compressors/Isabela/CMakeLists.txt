option(CBENCH_ENABLE_ISABELA "Enable ISABELA." ON)

if (CBENCH_ENABLE_ISABELA)
	set(MY_LIBS "")

    #Asks for Isabela library + set that in EXTRA_LIBS
    find_library(ISABELA_LIBRARY ${ISABELA_LIBRARY} CACHE STRING "ISABELA library")
    if (NOT ISABELA_LIBRARY)
    	message(FATAL_ERROR "Could NOT find ISABELA")
    else()
    	set(MY_LIBS ${MY_LIBS} ${ISABELA_LIBRARY})

		#Find other libraries
		#find_library(GSL_LIBRARY ${GSL_LIBRARY} CACHE STRING "GSL library")
		find_package(GSL REQUIRED)
    	set(MY_LIBS ${MY_LIBS} ${GSL_LIBRARIES})

    	#Asks for ISABELA include files + includes it in the project
    	set (ISABELA_INCLUDE_PATH ${ISABELA_INCLUDE_PATH} CACHE PATH "ISABELA include directory")

    	set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_ISABELA")
    	set_property(GLOBAL APPEND PROPERTY includeDirectories ${ISABELA_INCLUDE_PATH})
		set_property(GLOBAL APPEND PROPERTY includeDirectories ${GSL_INCLUDE_DIRS})
    	set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()