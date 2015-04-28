cmake_minimum_required(VERSION 2.8)

project(228a_homework)

SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Weffc++ -O3 -funroll-all-loops") 

# Path to hdf5 library 
find_library(hdf5 libhdf5.so.7 PATHS /usr/lib/)
find_library(hdf5_cpp libhdf5_cpp.so.7 PATHS /usr/lib/)

# Path to installation directory
set(MY_DIRECTORY .)

include_directories(${MY_DIRECTORY}/include)

set(SOURCES main.cpp)

add_executable(compute_task ${SOURCES})

target_link_libraries(compute_task ${hdf5} ${hdf5_cpp})

install(TARGETS compute_task RUNTIME DESTINATION bin)

set(CMAKE_BUILD_TYPE Release)
