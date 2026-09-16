cmake_minimum_required(VERSION 3.18)

# Compiler selection and macOS targeting must precede project().
option(DTCC_USE_HOMEBREW_LLVM "Use the Homebrew LLVM toolchain on macOS" OFF)
if(DTCC_USE_HOMEBREW_LLVM)
  if(NOT CMAKE_HOST_APPLE)
    message(FATAL_ERROR "DTCC_USE_HOMEBREW_LLVM requires macOS")
  endif()
  find_program(HOMEBREW_CLANGXX clang++
    PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin NO_DEFAULT_PATH)
  if(NOT HOMEBREW_CLANGXX)
    message(FATAL_ERROR "DTCC_USE_HOMEBREW_LLVM=ON but Homebrew clang++ was not found")
  endif()
  set(CMAKE_CXX_COMPILER "${HOMEBREW_CLANGXX}")
  message(STATUS "Using Homebrew LLVM toolchain (${HOMEBREW_CLANGXX})")
endif()

if(CMAKE_HOST_APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET
    AND NOT DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
  set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING "Minimum macOS version to target")
endif()

project(dtcc-core LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
endif()

option(DTCC_ENABLE_WARNINGS "Enable recommended warning flags" ON)
option(DTCC_USE_TRIANGLE "Enable the optional Triangle backend" OFF)
set(DTCC_TRIANGLE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dtcc_core/cpp/external/triangle"
  CACHE PATH "Directory or prefix containing the Triangle implementation header")

# Python and pybind11 are supplied by the isolated build environment.
set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_FRAMEWORK LAST)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE PYBIND11_DISCOVERY_RESULT
)
if(NOT PYBIND11_DISCOVERY_RESULT EQUAL 0)
  message(FATAL_ERROR "pybind11 must be installed in the Python build environment")
endif()
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
find_package(pybind11 CONFIG REQUIRED)

if(DTCC_USE_TRIANGLE)
  # Respect a changed prefix when reusing an incremental build directory.
  unset(TRIANGLE_INCLUDE_DIR CACHE)
  find_path(TRIANGLE_INCLUDE_DIR triangle.h
    HINTS ${DTCC_TRIANGLE_DIR} $ENV{TRIANGLE_DIR}
    PATH_SUFFIXES include include/triangle triangle
    NO_DEFAULT_PATH
  )
  if(NOT TRIANGLE_INCLUDE_DIR)
    message(FATAL_ERROR "DTCC_USE_TRIANGLE=ON but triangle.h was not found. Set DTCC_TRIANGLE_DIR to its directory or prefix.")
  endif()
  message(STATUS "Triangle enabled: include at ${TRIANGLE_INCLUDE_DIR}")
endif()

add_subdirectory(dtcc_core)
