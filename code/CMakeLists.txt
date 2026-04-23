cmake_minimum_required(VERSION 3.17)

project(TaylorSeriesComputeLog LANGUAGES C)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME} compute_log.c)