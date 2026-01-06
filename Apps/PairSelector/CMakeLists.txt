cmake_minimum_required(VERSION 3.16)
project(PairSelector VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(PairSelector main.cpp)

if(WIN32)
    target_compile_definitions(PairSelector PRIVATE _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
endif()

if(MSVC)
    target_compile_options(PairSelector PRIVATE /O2 /utf-8)
else()
    target_compile_options(PairSelector PRIVATE -O2)
endif()
