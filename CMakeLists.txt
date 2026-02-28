cmake_minimum_required(VERSION 3.12)
project(Lightcone LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(LIGHTCONE_COMMON_COMPILE_OPTIONS
  -march=native
  -mtune=native
  -flto
  -fno-stack-protector
  -Wall
  -Wextra
  -pedantic
)

add_executable(Lightcone Lightcone.cpp)
target_include_directories(Lightcone PRIVATE toml/include)
target_compile_options(Lightcone PRIVATE ${LIGHTCONE_COMMON_COMPILE_OPTIONS})
