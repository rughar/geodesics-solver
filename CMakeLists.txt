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

# Optional warnings for GCC / Clang
include(cmake/gcc-warnings.cmake OPTIONAL)

# ------------------------------------------------------------
# Tutorial / example binary
# ------------------------------------------------------------
add_executable(example_schwarzschild examples/schwarzschild_plane.cpp)
target_link_libraries(example_schwarzschild PRIVATE geodesics_solver)

# ------------------------------------------------------------
# Install rules (so users can do find_package(geodesics_solver))
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