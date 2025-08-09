cmake_minimum_required(VERSION 3.20)
project(geodesics-solver VERSION 0.1.0 LANGUAGES CXX)

# ------------------------------------------------------------
# Header-only library target
# ------------------------------------------------------------
add_library(geodesics_solver INTERFACE)
target_include_directories(geodesics_solver
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(geodesics_solver INTERFACE cxx_std_17)

# Optional warnings for GCC / Clang (if available)
include(cmake/gcc-warnings.cmake OPTIONAL)

# ------------------------------------------------------------
# Example binary
# ------------------------------------------------------------
add_executable(example_schwarzschild examples/schwarzschild_plane.cpp)
target_link_libraries(example_schwarzschild PRIVATE geodesics_solver)

# ------------------------------------------------------------
# Profiling-friendly build flags
#   - MSVC (Visual Studio): /O2 /Zi /Zo /Ob1 /Oy-  + linker /DEBUG:FULL /INCREMENTAL:NO
#   - GCC/Clang: -O2 -g -fno-omit-frame-pointer
# ------------------------------------------------------------
if(MSVC)
  # Compiler options for RelWithDebInfo
  target_compile_options(example_schwarzschild PRIVATE
    $<$<CONFIG:RelWithDebInfo>:/O2>
    $<$<CONFIG:RelWithDebInfo>:/Zi>
    $<$<CONFIG:RelWithDebInfo>:/Zo>     # enhanced optimized debug info
    $<$<CONFIG:RelWithDebInfo>:/Ob1>    # limit inlining
    $<$<CONFIG:RelWithDebInfo>:/Oy->    # keep frame pointers
  )
  # Linker options for RelWithDebInfo
  target_link_options(example_schwarzschild PRIVATE
    $<$<CONFIG:RelWithDebInfo>:/DEBUG:FULL>
    $<$<CONFIG:RelWithDebInfo>:/INCREMENTAL:NO>
  )
else()
  # GCC/Clang
  target_compile_options(example_schwarzschild PRIVATE
    $<$<CONFIG:RelWithDebInfo>:-O2>
    $<$<CONFIG:RelWithDebInfo>:-g>
    $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>
  )
endif()


# Set default build type to RelWithDebInfo for single-config generators
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# ------------------------------------------------------------
# Install rules
# ------------------------------------------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/geodesics_solverConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS geodesics_solver EXPORT geodesics_solverTargets)
install(EXPORT  geodesics_solverTargets
        FILE    geodesics_solverTargets.cmake
        NAMESPACE geodesics::
        DESTINATION lib/cmake/geodesics_solver)
install(FILES "${PROJECT_BINARY_DIR}/geodesics_solverConfigVersion.cmake"
        DESTINATION lib/cmake/geodesics_solver)