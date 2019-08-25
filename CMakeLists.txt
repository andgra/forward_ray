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

add_executable(forward_ray src/main.cpp src/point.h src/vec2.h src/data.h src/edge.h src/beam.h src/intersection.h src/figure.h src/arrayOfReceiversTransmitters.h src/algo.h src/dft.h src/comp.cpp src/lib/fft.cpp src/helper.h src/exo.h src/speedmap.h src/lib/bitmap_image.hpp)

target_link_libraries(forward_ray ${MPI_LIBRARY})

set(INCLUDE_MPI_DLL false)

if (INCLUDE_MPI_DLL)
    #msmpi.dll
    add_custom_command(TARGET forward_ray POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/src/lib/msmpi.dll      # <--this is in-file
            $<TARGET_FILE_DIR:forward_ray>)                 # <--this is out-file path
endif()