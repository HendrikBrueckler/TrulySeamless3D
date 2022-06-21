set(CLI_NAME "TrulySeamless3D")

add_executable(${CLI_NAME} main.cc)
target_link_libraries(${CLI_NAME} TS3D::TS3D)
set_target_properties(${CLI_NAME} PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON)
