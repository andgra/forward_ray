cmake_minimum_required(VERSION 3.13)
project(straight_ray)

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

configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")


set(FFTW_ROOT "D:\\projects\\clion\\fftw-3.3.5-dll32")
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
link_directories (${FFTW_INCLUDE_DIRS})
#message(${FFTW_LIBRARIES})


#set(OpenCV_DIR "D:\\opencv\\build\\")
#set(OpenCV_DIR "D:\\opencv\\build\\x64\\vc14\\lib")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
#include_directories(${OpenCV_INCLUDE_DIRS})
#find_package(OpenCV REQUIRED)
#include_directories(D:\\opencv\\build\\)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static -lpthread")

add_executable(straight_ray main.cpp point.h vec2.h data.h edge.h beam.h intersection.h figure.h arrayOfReceiversTransmitters.h algo.h dft.h comp.cpp lib/fft.cpp helper.h exo.h speedmap.h lib/bitmap_image.hpp)

target_link_libraries(straight_ray ${FFTW_LIBRARIES})

add_custom_command(TARGET straight_ray POST_BUILD        # Adds a post-build event to MyTest
        COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
        ${FFTW_LIBRARIES}      # <--this is in-file
        $<TARGET_FILE_DIR:straight_ray>)                 # <--this is out-file path
#message(${MPI_LIBRARY})
#target_link_libraries(straight_ray ${MPI_LIBRARY} "-Wl,--enable-stdcall-fixup")
#target_link_libraries(straight_ray ${MPI_LIBRARY} "-static_mpi")
target_link_libraries(straight_ray ${MPI_LIBRARY})

#msmpi.dll
add_custom_command(TARGET straight_ray POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/lib/msmpi.dll      # <--this is in-file
        $<TARGET_FILE_DIR:straight_ray>)                 # <--this is out-file path
