cmake_minimum_required(VERSION 3.18)

project(dtcc-core LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
  message(STATUS "Build type not specified, defaulting to Release")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(DTCC_ENABLE_WARNINGS "Enable recommended warning flags" ON)
option(DTCC_USE_SPADE "Enable SPADE triangulation" ON)
option(DTCC_USE_HOMEBREW_LLVM "Use the Homebrew LLVM toolchain on macOS" OFF)

if(APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING "Minimum macOS version to target" FORCE)
endif()

if(APPLE AND DTCC_USE_HOMEBREW_LLVM)
  find_program(HOMEBREW_CLANG /opt/homebrew/opt/llvm/bin/clang)
  find_program(HOMEBREW_CLANGXX /opt/homebrew/opt/llvm/bin/clang++)
  if(HOMEBREW_CLANG AND HOMEBREW_CLANGXX)
    set(CMAKE_C_COMPILER ${HOMEBREW_CLANG})
    set(CMAKE_CXX_COMPILER ${HOMEBREW_CLANGXX})
    message(STATUS "Using Homebrew LLVM toolchain (${HOMEBREW_CLANGXX})")
    find_path(HOMEBREW_LIBOMP_PREFIX include/omp.h
      PATHS /opt/homebrew/opt/libomp /usr/local/opt/libomp
    )
    if(HOMEBREW_LIBOMP_PREFIX)
      list(APPEND CMAKE_PREFIX_PATH ${HOMEBREW_LIBOMP_PREFIX})
    endif()
  else()
    message(FATAL_ERROR "DTCC_USE_HOMEBREW_LLVM=ON but clang++ not found in /opt/homebrew/opt/llvm")
  endif()
endif()

if(DTCC_ENABLE_WARNINGS)
  if(MSVC)
    add_compile_options(/W4 /permissive-)
  else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-reorder)
  endif()
endif()

if(NOT MSVC)
  add_compile_options(-fvisibility=hidden)
endif()

# ------------------------------------------------------------------------------
# Python and pybind11 discovery
# ------------------------------------------------------------------------------
set(Python_FIND_VIRTUALENV FIRST)
set(Python_FIND_FRAMEWORK LAST)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(PYBIND11_CMAKE_DIR)
  list(APPEND CMAKE_PREFIX_PATH ${PYBIND11_CMAKE_DIR})
endif()
find_package(pybind11 CONFIG REQUIRED)

# ------------------------------------------------------------------------------
# Include directories and dependencies
# ------------------------------------------------------------------------------
set(DTCC_CPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dtcc_core/cpp)
set(DTCC_CORE_PUBLIC_INCLUDE_DIRS
  ${DTCC_CPP_DIR}/include
  ${DTCC_CPP_DIR}/include/citymodel
)
set(DTCC_CORE_SYSTEM_INCLUDE_DIRS
  ${DTCC_CPP_DIR}/external
  ${DTCC_CPP_DIR}/external/triangle
  ${DTCC_CPP_DIR}/external/amgcl
)

find_package(Boost)
if(Boost_FOUND)
  message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
  list(APPEND DTCC_CORE_SYSTEM_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
  set(DTCC_HAVE_BOOST TRUE)
else()
  message(STATUS "Boost not found; building AMGCL without Boost support")
  set(DTCC_HAVE_BOOST FALSE)
endif()

# ------------------------------------------------------------------------------
# SPADE integration (sibling checkout)
# ------------------------------------------------------------------------------
if(DTCC_USE_SPADE)
  set(DTCC_SPADE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../dtcc-wrapper-spade/cppspade"
    CACHE PATH "Location of the dtcc-wrapper-spade/cppspade directory"
  )

  if(EXISTS "${DTCC_SPADE_DIR}/CMakeLists.txt")
    find_program(CARGO_EXECUTABLE cargo)
    if(CARGO_EXECUTABLE)
      set(DTCC_SPADE_INCLUDE_DIR "${DTCC_SPADE_DIR}/include"
        CACHE PATH "SPADE wrapper include directory")
      set(DTCC_SPADE_WRAPPER_SRC "${DTCC_SPADE_DIR}/src/spade_wrapper.cpp"
        CACHE FILEPATH "SPADE wrapper source file")

      if(APPLE)
        set(DTCC_SPADE_STATIC_LIB "${DTCC_SPADE_DIR}/target/release/libspade_ffi.a")
      elseif(UNIX)
        set(DTCC_SPADE_STATIC_LIB "${DTCC_SPADE_DIR}/target/release/libspade_ffi.a")
      else()
        set(DTCC_SPADE_STATIC_LIB "${DTCC_SPADE_DIR}/target/release/spade_ffi.lib")
      endif()
      set(DTCC_SPADE_STATIC_LIB ${DTCC_SPADE_STATIC_LIB}
        CACHE FILEPATH "Prebuilt SPADE static library")

      add_custom_command(
        OUTPUT ${DTCC_SPADE_STATIC_LIB}
        COMMAND ${CARGO_EXECUTABLE} build --release
        WORKING_DIRECTORY ${DTCC_SPADE_DIR}
        COMMENT "Building Rust SPADE FFI library"
        VERBATIM
      )
      add_custom_target(spade_ffi_build DEPENDS ${DTCC_SPADE_STATIC_LIB})
      message(STATUS "SPADE integration enabled (${DTCC_SPADE_DIR})")
    else()
      message(WARNING "cargo not available; disabling SPADE integration")
      set(DTCC_USE_SPADE OFF CACHE BOOL "" FORCE)
    endif()
  else()
    message(WARNING "dtcc-wrapper-spade not detected at ${DTCC_SPADE_DIR}; disabling SPADE integration")
    set(DTCC_USE_SPADE OFF CACHE BOOL "" FORCE)
  endif()
endif()

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY ON)

add_subdirectory(dtcc_core)
