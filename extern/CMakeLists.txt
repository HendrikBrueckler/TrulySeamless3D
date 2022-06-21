find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    git_update_submodules("${CMAKE_CURRENT_SOURCE_DIR}")
else()
    message(WARNING "TRULYSEAMLESS_SUBMODULES_MANUAL is on, but this is not a git repository")
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/libHexEx/CMakeLists.txt")
    if (${TRULYSEAMLESS_STANDALONE})
        message(FATAL_ERROR "The submodules were not downloaded! Please update submodules and try again.")
    else()
        message(WARNING "TrulySeamless3D submodules not downloaded! Must supply them in a base project.")
    endif()
endif()
