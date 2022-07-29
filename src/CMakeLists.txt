### Find dependencies
# OVM
if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh extern/OpenVolumeMesh EXCLUDE_FROM_ALL)
    else()
        find_package(OpenVolumeMesh REQUIRED)
    endif()
endif()
list(APPEND TS3D_LIB_LIST OpenVolumeMesh::OpenVolumeMesh)

# libhexex
if (NOT TARGET HexEx::HexEx)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/libHexEx/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/libHexEx extern/libHexEx EXCLUDE_FROM_ALL)
    else()
        find_package(HexEx REQUIRED)
    endif()
endif()
list(APPEND TS3D_LIB_LIST HexEx::HexEx)

#eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND TS3D_LIB_LIST Eigen3::Eigen)

### Source files
list(APPEND TS3D_SOURCE_LIST
    "helpers.cc"
    "timer.cc"
    "trulyseamless.cc")
list(APPEND TS3D_HEADER_LIST
    "../include/TS3D/helpers.h"
    "../include/TS3D/timer.h"
    "../include/TS3D/trulyseamless.h")

### Create target
add_library(TS3D ${TS3D_SOURCE_LIST} ${TS3D_HEADER_LIST})
add_library(TS3D::TS3D ALIAS TS3D)

### Configure install paths
include(CMakePackageConfigHelpers)
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR "CMake")
else()
  set(DEF_INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/TS3D")
endif()
set(TS3D_INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH "Installation directory for CMake files")
set(TS3D_INSTALL_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Installation directory for libraries")
set(TS3D_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Installation directory for header files")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../include")
set(CONF_LIB_DIRS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
configure_package_config_file(../TS3DConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/TS3DConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
# ... for the install tree
set(CONF_INCLUDE_DIRS "${TS3D_INSTALL_INCLUDE_DIR}")
set(CONF_LIB_DIRS "${TS3D_INSTALL_LIB_DIR}")
configure_package_config_file(../TS3DConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TS3DConfig.cmake"
  INSTALL_DESTINATION "${TS3D_INSTALL_CMAKE_DIR}")

### Include own headers with public access
target_include_directories(TS3D
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>"
                           "$<INSTALL_INTERFACE:${TS3D_INSTALL_INCLUDE_DIR}>")

### Link with dependencies
target_link_libraries_system(TS3D PUBLIC ${TS3D_LIB_LIST})

### Set properties
# properties
set_target_properties(TS3D PROPERTIES
                      DEBUG_POSTFIX "d"
                      VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      PUBLIC_HEADER "${TS3D_HEADER_LIST}"
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON)
# compile options
list(APPEND TS3D_COMPILE_OPTIONS_PRV "-Wall" "-Wextra" "-Wshadow" "-Wunused" "-Wunused-function")
target_compile_options(TS3D PRIVATE ${TS3D_COMPILE_OPTIONS_PRV})
# preprocessor defines
if (NOT ${TS3D_ENABLE_LOGGING})
    list(APPEND TS3D_COMPILE_DEFINITIONS_PRV "TRULYSEAMLESS_SILENT")
endif()
target_compile_definitions(TS3D PRIVATE ${TS3D_COMPILE_DEFINITIONS_PRV})

### Install and export
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TS3DConfig.cmake"
        DESTINATION "${TS3D_INSTALL_CMAKE_DIR}")

install(TARGETS       TS3D
        EXPORT        TS3DTargets
        LIBRARY       DESTINATION "${TS3D_INSTALL_LIB_DIR}"
        ARCHIVE       DESTINATION "${TS3D_INSTALL_LIB_DIR}"
        PUBLIC_HEADER DESTINATION "${TS3D_INSTALL_INCLUDE_DIR}/TS3D")

install(EXPORT TS3DTargets
        FILE TS3DTargets.cmake
        NAMESPACE TS3D::
        DESTINATION "${TS3D_INSTALL_CMAKE_DIR}")

export(EXPORT TS3DTargets
       NAMESPACE TS3D::)
