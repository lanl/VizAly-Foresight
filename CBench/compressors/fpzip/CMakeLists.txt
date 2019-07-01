option(CBENCH_ENABLE_FPZIP "Enable FPZIP." OFF)

if (CBENCH_ENABLE_FPZIP)
    set(MY_LIBS "")

    #Asks for fpzip library + set that in EXTRA_LIBS
    find_library(FPZIP_LIBRARY ${FPZIP_LIBRARY} CACHE STRING "FPZIP library")
    if (NOT FPZIP_LIBRARY)
        message(FATAL_ERROR "Could NOT find FPZIP")
    else()
        set(MY_LIBS ${MY_LIBS} ${FPZIP_LIBRARY})

        #Asks for fpzip include files + includes it in the project
        set (FPZIP_INCLUDE_PATH ${FPZIP_INCLUDE_PATH} CACHE PATH "FPZIP include directory")

        set_property(GLOBAL APPEND PROPERTY definitions "-DCBENCH_HAS_FPZIP")
        set_property(GLOBAL APPEND PROPERTY includeDirectories ${FPZIP_INCLUDE_PATH})
        set_property(GLOBAL APPEND PROPERTY linkLibraries ${MY_LIBS})
    endif()
endif()