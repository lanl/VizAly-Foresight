# library CMakeLists.txt template

# project name
project(genericio)

# captial project name string
string(TOUPPER ${PROJECT_NAME} NAME_ID)

# option for STATIC or SHARED
option(${NAME_ID}_BUILD_STATIC
	"Build a static version of this library"
	OFF
)

# project specific options
option(${NAME_ID}_ENABLE_BGQ_ENVIRONMENT
	"Enable BGQ Environment"
	OFF
)

# add source files
set(PROJECT_SOURCE_FILES
	DataFilter.cxx
	GenericIO.cxx
	GenericIOBase.cxx
	GenericIOMPIReader.cxx
	GenericIOPosixReader.cxx
	GenericIOReader.cxx
	GenericIOUtilities.cxx
)

# add header files
set(PROJECT_HEADER_FILES
#	include/source.h
)

# add interface header files
set(PROJECT_INTERFACE_HEADER_FILES
	include/AssocVector.h
	include/CRC64.h
	include/DataFilter.h
	include/GenericIOBase.h
	include/GenericIOCompress.h
	include/GenericIODefinitions.hpp
	include/GenericIO.h
	include/GenericIOMPIReader.h
	include/GenericIOPosixReader.h
	include/GenericIOReader.h
	include/GenericIOUtilities.h
)

# add library targets
if(${NAME_ID}_BUILD_STATIC)
	add_library(${PROJECT_NAME} STATIC
		${PROJECT_SOURCE_FILES}
		${PROJECT_HEADER_FILES}
		${PROJECT_INTERFACE_HEADER_FILES}
	)
else()
	add_library(${PROJECT_NAME} SHARED
		${PROJECT_SOURCE_FILES}
		${PROJECT_HEADER_FILES}
		${PROJECT_INTERFACE_HEADER_FILES}
	)
endif()

# gather compile definitions
set(DEFINITIONS __STDC_CONSTANT_MACROS)

if(${NAME_ID}_ENABLE_BGQ_ENVIRONMENT)
	set(DEFINITIONS ${DEFINITIONS} __bgq__)
endif()

# add compile definitions
# if no defintions need to be set, this command needs to be commented out
target_compile_definitions(${PROJECT_NAME}
	PUBLIC ${DEFINITIONS}
)

# add compile flags
# if no definitions need to be set, this command needs to be commented out
#target_compile_options(${PROJECT}
#)

# add include directories
target_include_directories(${PROJECT_NAME}
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# add library dependencies
target_link_libraries(${PROJECT_NAME}
	PUBLIC mpi
	PUBLIC openmp
)

# add doxygen API documentation
find_package(Doxygen)

if(DOXYGEN_FOUND)
	configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		@ONLY
	)

	add_custom_target(doxygen
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen"
		VERBATIM
	)
endif()

# install instructions
if(${NAME_ID}_BUILD_STATIC)
	install(TARGETS ${PROJECT_NAME}
		DESTINATION lib
	)
else()
	if(WIN32)
		install(TARGETS ${PROJECT_NAME}
			DESTINATION bin
		)
	else()
		install(TARGETS ${PROJECT_NAME}
			DESTINATION lib/3rdparty/${PROJECT_NAME}
		)
	endif()
endif()

install(FILES
	${PROJECT_INTERFACE_HEADER_FILES}
	DESTINATION include/${PROJECT_NAME}
)

#install(TARGETS doxygen
#	DESTINATION docs
#)
