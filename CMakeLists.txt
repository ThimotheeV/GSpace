cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_CXX_STANDARD 17)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

project(simulateur VERSION 0.0.1 LANGUAGES CXX)

include(CTest)
include(lib/Catch/Catch.cmake)

enable_testing()

option(CMAKE_FUNC_TEST OFF)
# pour pouvoir lui donner des options de compilation (e.g. dans sources/CMakeLists.txt : target_compile_options() )
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

#add path search for #include <> in CXX files for the compiler
include_directories ("${PROJECT_SOURCE_DIR}/lib/Glib")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Discret_Sichel")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Catch")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Tskit")
include_directories ("${PROJECT_SOURCE_DIR}/sources")
include_directories(BEFORE SYSTEM "${PROJECT_SOURCE_DIR}")

#include subprojet in main cmake projet, so that cmake can find other cmake files.

add_subdirectory(lib/Glib)
add_subdirectory(lib/Discret_Sichel)
add_subdirectory(lib/Tskit)
add_subdirectory(sources)

if (CMAKE_FUNC_TEST)
    add_subdirectory(func_tests)
else()
    add_subdirectory(tests)
endif()

# fait le .o
add_executable(GSpace main.cpp)
# lie les .o/.a ...
target_link_libraries(GSpace libsimu disSichel libcalc tskit)
