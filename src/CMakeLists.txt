include(macrOS)

### Better safe than sorry...
list(APPEND ${TRULYSEAMLESS_LIB_NAME}_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror" "-Wshadow" "-Wunused" "-Wunused-function")

### Compile time options
if (NOT ${TRULYSEAMLESS_ENABLE_LOGGING})
    list(APPEND ${TRULYSEAMLESS_LIB_NAME}_COMPILE_DEFINITIONS "TRULYSEAMLESS_SILENT")
endif()

### Add all source files
list(APPEND ${TRULYSEAMLESS_LIB_NAME}_SOURCE_LIST
    "helpers.cc"
    "timer.cc"
    "trulyseamless.cc"
)
add_library(${TRULYSEAMLESS_LIB_NAME} STATIC ${${TRULYSEAMLESS_LIB_NAME}_SOURCE_LIST})

### Find dependencies
# OVM
if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh extern/OpenVolumeMesh EXCLUDE_FROM_ALL)
endif()
list(APPEND ${TRULYSEAMLESS_LIB_NAME}_LIB_LIST OpenVolumeMesh::OpenVolumeMesh)

# libhexex
if (NOT TARGET HexEx::HexEx)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/libHexEx extern/libHexEx EXCLUDE_FROM_ALL)
endif()
list(APPEND ${TRULYSEAMLESS_LIB_NAME}_LIB_LIST HexEx::HexEx)

#eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND ${TRULYSEAMLESS_LIB_NAME}_LIB_LIST Eigen3::Eigen)

### Include own headers with public access
target_include_directories(${TRULYSEAMLESS_LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

### Include headers as system headers to prevent warnings/errors from strict warnings
target_include_directories(${TRULYSEAMLESS_LIB_NAME} SYSTEM PUBLIC ${${TRULYSEAMLESS_LIB_NAME}_PUB_INC_LIST})

### Link with dependencies
target_link_libraries_system(${TRULYSEAMLESS_LIB_NAME} PUBLIC ${${TRULYSEAMLESS_LIB_NAME}_LIB_LIST})

### Set compile options
set_target_properties(${TRULYSEAMLESS_LIB_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_options(${TRULYSEAMLESS_LIB_NAME} PRIVATE ${${TRULYSEAMLESS_LIB_NAME}_COMPILE_OPTIONS})
target_compile_definitions(${TRULYSEAMLESS_LIB_NAME} PRIVATE ${${TRULYSEAMLESS_LIB_NAME}_COMPILE_DEFINITIONS})
