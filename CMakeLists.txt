cmake_minimum_required(VERSION 3.13)
project(forward_ray)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set(MPI_C_HEADER_DIR "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include")
set(MPI_CXX_HEADER_DIR "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include")
set(MPI_CXX_FLAGS "-static_mpi")
include(FindMPI)
FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static -lpthread")

add_executable(forward_ray main.cpp point.h vec2.h data.h edge.h beam.h intersection.h figure.h arrayOfReceiversTransmitters.h algo.h dft.h comp.cpp lib/fft.cpp helper.h exo.h speedmap.h lib/bitmap_image.hpp)

target_link_libraries(forward_ray ${MPI_LIBRARY})

#msmpi.dll
add_custom_command(TARGET forward_ray POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/lib/msmpi.dll      # <--this is in-file
        $<TARGET_FILE_DIR:forward_ray>)                 # <--this is out-file path
