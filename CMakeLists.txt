cmake_minimum_required(VERSION 3.12)
project(LightCone CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-march=native -mtune=native -flto -fno-stack-protector -Wall -Wextra -pedantic)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-stack_size -Wl,0x1000000")
endif()

# Add toml include directory
include_directories(toml/include)

# Add executable
add_executable(LightCone LightCone.cpp)
