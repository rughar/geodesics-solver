if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(geodesics_solver INTERFACE
      -Wall -Wextra -Wpedantic -Wconversion -Wno-sign-conversion)
endif()