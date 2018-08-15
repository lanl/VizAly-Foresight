# root CMakeLists.txt file

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6.2)

# solution name
project(case_study)

# set cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/")

# 3rdparty dependencies
#add_subdirectory(3rdparty/c-blosc-1.10.2)
#add_subdirectory(3rdparty/fpzip-1.1.0)
add_subdirectory(3rdparty/genericio)
#add_subdirectory(3rdparty/sz-1.4.9-hacc)
#add_subdirectory(3rdparty/z-checker-installer/SZ)
#add_subdirectory(3rdparty/xz-5.2.3)
#add_subdirectory(3rdparty/zfp-0.5.0)
add_subdirectory(3rdparty/mpi)
#add_subdirectory(3rdparty/openmp)
add_subdirectory(3rdparty/BigCrunch)

# utility project
add_subdirectory(utility)

# data loader project
add_subdirectory(data_loader)

# benchmark projects
add_subdirectory(benchmarks/blosc)
#add_subdirectory(benchmarks/gio)
#add_subdirectory(benchmarks/lzma)
#add_subdirectory(benchmarks/fpzip)
#add_subdirectory(benchmarks/sz)
#add_subdirectory(benchmarks/zfp)
add_subdirectory(benchmarks/BigCrunch)

# high-level install instructions (e.g. scripts, resources, data, ...)
install(DIRECTORY
	${CMAKE_CURRENT_SOURCE_DIR}/../env/
	DESTINATION env
)

# install run_all.sh script
install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/../run/run_all.sh
	DESTINATION run
)

# create folders
install(DIRECTORY
	DESTINATION run/logs
)
