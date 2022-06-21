cmake_minimum_required(VERSION 3.1)

project(TS3D VERSION 1.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(macrOS)
include(GNUInstallDirs)

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

### Standalone?
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(TS3D_STANDALONE TRUE)
    message("Configuring TS3D as standalone project...")
else()
    set(TS3D_STANDALONE FALSE)
    message("Configuring TS3D inside another cmake project...")
endif()

### Options
option(TS3D_BUILD_CLI         "Build CLI app for TrulySeamless3D" ${TS3D_STANDALONE})
option(TS3D_ENABLE_LOGGING    "Allow logging of TrulySeamless3D"  ${TS3D_STANDALONE})
option(TS3D_SUBMODULES_MANUAL "Skip automatic submodule download" OFF)
option(BUILD_SHARED_LIBS      "Build libraries as shared as opposed to static" ON)

### Global cmake settings
set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_EXTENSIONS OFF)
if (TS3D_STANDALONE)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
endif()

### Load external libs
if(NOT ${TS3D_SUBMODULES_MANUAL})
    add_subdirectory(extern)
endif()

### Add own lib
# lib target: TS3D::TS3D
add_subdirectory(src)

### Add CLI
if (TS3D_BUILD_CLI)
    add_subdirectory(cli)
endif()

### Fake successful finder run if compiling as a dependent project.
if (NOT TS3D_STANDALONE)
    set(TS3D_FOUND true PARENT_SCOPE)
    set(TS3D_LIBRARIES TS3D PARENT_SCOPE)
    set(TS3D_LIBRARY TS3D PARENT_SCOPE)
    set(TS3D_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include/" PARENT_SCOPE)
    set(TS3D_LIBRARY_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" CACHE PATH "The directory where the TS3D libraries can be found.")
endif()
