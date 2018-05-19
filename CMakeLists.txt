cmake_minimum_required(VERSION 3.10)
project(straightRay)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static -lpthread")

add_executable(straightRay main.cpp)