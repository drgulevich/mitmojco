cmake_minimum_required(VERSION 3.5.1)
project(mitmojco VERSION 1.1)
set(CMAKE_BUILD_TYPE Release)

include_directories(include)

set(SOURCES src/mitmojco.c src/opt_filter.c)

add_library(mitmojco SHARED ${SOURCES})

set_target_properties(mitmojco PROPERTIES 
    VERSION ${PROJECT_VERSION}
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
