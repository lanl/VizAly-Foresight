
if (CBENCH_ENABLE_FPZIP)

    find_library(FPZIP_LIBRARY ${FPZIP_LIBRARY} CACHE STRING "FPZIP library")

    if (FPZIP_LIBRARY)
        set(FPZIP_INCLUDE_PATH ${FPZIP_INCLUDE_PATH} CACHE PATH "FPZIP include directory")
		target_compile_definitions(CBench PRIVATE -DCBENCH_HAS_FPZIP)
		target_include_directories(CBench PRIVATE ${FPZIP_INCLUDE_PATH})
        target_include_directories(GenericIO PRIVATE ${FPZIP_INCLUDE_PATH})
        set_property(TARGET CBench APPEND PROPERTY LIB "${FPZIP_LIBRARY}")
    else()
		message(FATAL_ERROR "Could NOT find FPZIP")
    endif()
endif()