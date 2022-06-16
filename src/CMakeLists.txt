include(macrOS)

### Better safe than sorry...
list(APPEND TS3D_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror" "-Wshadow" "-Wunused" "-Wunused-function")

### Compile time options
if (NOT ${TS3D_ENABLE_LOGGING})
    list(APPEND TS3D_COMPILE_DEFINITIONS "TRULYSEAMLESS_SILENT")
endif()

### Add all source files
list(APPEND TS3D_SOURCE_LIST
    "helpers.cc"
    "timer.cc"
    "trulyseamless.cc"
)
add_library(${TS3D_LIB_NAME} STATIC ${TS3D_SOURCE_LIST})

### Find dependencies
# OVM
if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh extern/OpenVolumeMesh EXCLUDE_FROM_ALL)
endif()
list(APPEND TS3D_LIB_LIST OpenVolumeMesh::OpenVolumeMesh)

# libhexex
if (NOT TARGET HexEx::HexEx)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/libHexEx extern/libHexEx EXCLUDE_FROM_ALL)
endif()
list(APPEND TS3D_LIB_LIST HexEx::HexEx)

#eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND TS3D_LIB_LIST Eigen3::Eigen)

### Include own headers with public access
target_include_directories(${TS3D_LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

### Include headers as system headers to prevent warnings/errors from strict warnings
target_include_directories(${TS3D_LIB_NAME} SYSTEM PUBLIC ${TS3D_PUB_INC_LIST})

### Link with dependencies
target_link_libraries_system(${TS3D_LIB_NAME} PUBLIC ${TS3D_LIB_LIST})

### Set compile options
set_target_properties(${TS3D_LIB_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_options(${TS3D_LIB_NAME} PRIVATE ${TS3D_COMPILE_OPTIONS})
target_compile_definitions(${TS3D_LIB_NAME} PRIVATE ${TS3D_COMPILE_DEFINITIONS})
