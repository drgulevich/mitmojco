cmake_minimum_required(VERSION 2.6)
project(mitmojco)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS_RELEASE "-Ofast -march=native -fopenmp")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall -fopenmp")

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
if(${CMAKE_BUILD_TYPE} STREQUAL Release)
    message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL Debug)
    message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
endif()

include_directories(include)

set(SOURCES src/mitmojco.c src/opt_filter.c)

add_library(mitmojco SHARED ${SOURCES})

set_target_properties(mitmojco PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    PUBLIC_HEADER "include/mitmojco/mitmojco.h;include/mitmojco/opt_filter.h")

install(TARGETS mitmojco
    LIBRARY DESTINATION /usr/local/lib
    PUBLIC_HEADER DESTINATION /usr/local/include/mitmojco)

add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
    COMMAND ${CMAKE_COMMAND} -E remove install_manifest.txt *.so.*
    COMMENT "distclean invoked"
    )
