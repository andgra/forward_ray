cmake_minimum_required(VERSION 3.10)
project(straight_ray)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set(MPI_C_HEADER_DIR "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include")
set(MPI_CXX_HEADER_DIR "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include")
include(FindMPI)
FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

#set(OpenCV_DIR "D:\\opencv\\build\\")
#set(OpenCV_DIR "D:\\opencv\\build\\x64\\vc14\\lib")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
#include_directories(${OpenCV_INCLUDE_DIRS})
#find_package(OpenCV REQUIRED)
#include_directories(D:\\opencv\\build\\)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static -lpthread")

add_executable(straight_ray main.cpp point.h vec2.h data.h edge.h beam.h intersection.h figure.h arrayOfReceiversTransmitters.h algo.h dft.h comp.cpp fft.cpp helper.h exo.h)

target_link_libraries(straight_ray ${MPI_LIBRARY} "-Wl,--enable-stdcall-fixup")
