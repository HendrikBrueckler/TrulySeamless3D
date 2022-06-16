cmake_minimum_required(VERSION 3.1)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(macrOS)

project(TrulySeamless3D VERSION 1.0)

### Default to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build this in a source directory (or any directory with a CMakeLists.txt file). "
                        "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

### Build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "RelWithDebInfo")
endif()

### Binary output
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/bin/)

### Standalone?
if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(TS3D_STANDALONE TRUE)
else()
    set(TS3D_STANDALONE FALSE)
endif()

### Options
option(TS3D_BUILD_CLI         "Build CLI app for TrulySeamless3D" ${TS3D_STANDALONE})
option(TS3D_ENABLE_LOGGING    "Allow logging of TrulySeamless3D"  ${TS3D_STANDALONE})
option(TS3D_SUBMODULES_MANUAL "Skip automatic submodule download" OFF)

### Load external libs
add_subdirectory(extern)

### Add own lib
set(TS3D_LIB_NAME trulyseamless3d)
add_subdirectory(src)

### Add CLI
if (TS3D_BUILD_CLI)
    add_subdirectory(cli)
endif()
