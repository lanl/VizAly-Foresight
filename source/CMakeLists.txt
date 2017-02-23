# root CMakeLists.txt file

# minimum cmake version requirement
cmake_minimum_required(VERSION 3.6.2)

# solution name
project(case_study)

# 3rdparty dependencies
add_subdirectory(3rdparty/c-blosc-1.10.2)
add_subdirectory(3rdparty/fpzip-1.1.0)
add_subdirectory(3rdparty/genericio)
#add_subdirectory(3rdparty/sz-1.4.8-beta)
add_subdirectory(3rdparty/sz-1.4.9-hacc)
add_subdirectory(3rdparty/xz-5.2.3)
add_subdirectory(3rdparty/zfp-0.5.0)
add_subdirectory(3rdparty/mpi)
add_subdirectory(3rdparty/openmp)
#add_subdirectory(3rdparty/netcdf-cxx4-4.3.0)

# utility project
add_subdirectory(utility)

# data loader project
add_subdirectory(data_loader)

# benchmark projects
add_subdirectory(benchmarks/blosc)
add_subdirectory(benchmarks/gio)
add_subdirectory(benchmarks/lzma)
add_subdirectory(benchmarks/fpzip)
add_subdirectory(benchmarks/sz)
add_subdirectory(benchmarks/zfp)

# install instructions
install(DIRECTORY
	${CMAKE_CURRENT_SOURCE_DIR}/../env/
	DESTINATION env
)

install(DIRECTORY
	DESTINATION run/logs
)
