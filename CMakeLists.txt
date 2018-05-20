cmake_minimum_required(VERSION 3.10)
project(straight_ray)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static -lpthread")

add_executable(straight_ray main.cpp point.h vec2.h edge.h #[[data.h beam.h intersection.h figure.h arrayOfReceiversTransmitters.h algo.h dft.h dft.cpp dft.h]] diffs.h)