# Project configuration
cmake_minimum_required(VERSION 3.15)
project(dtcc-core)

# Note: To building in Debug mode, use cmake .. -DCMAKE_BUILD_TYPE=Debug
# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type not specified, using Release by default")
endif()

option(DTCC_USE_HOMEBREW_LLVM "Use the Homebrew LLVM toolchain on macOS" OFF)

# MacOS doesn't ship with OpenMP, so we need to use Homebrew's Clang
# if we want OpenMP support. To get OpenMP support on Mac run:
# brew install llvm
# brew install libomp
if(DTCC_USE_HOMEBREW_LLVM AND APPLE)
  find_program(HOMEBREW_CLANG /opt/homebrew/opt/llvm/bin/clang)
  find_program(HOMEBREW_CLANGXX /opt/homebrew/opt/llvm/bin/clang++)

  # Set Homebrew compiler paths, if found, and set OpenMP flags
  if(HOMEBREW_CLANG AND HOMEBREW_CLANGXX)
    # set(CMAKE_PREFIX_PATH $(brew --prefix libomp))
    # Get Homebrew prefix path for libomp
    execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX)
    string(STRIP ${HOMEBREW_LIBOMP_PREFIX} HOMEBREW_LIBOMP_PREFIX)
    list(APPEND CMAKE_PREFIX_PATH ${HOMEBREW_LIBOMP_PREFIX})
    message(STATUS "Using Homebrew Clang (${HOMEBREW_CLANG}) and Clang++ (${HOMEBREW_CLANGXX}) compilers")
    message(STATUS "CMake prefix ${CMAKE_PREFIX_PATH}")
    set(CMAKE_C_COMPILER ${HOMEBREW_CLANG})
    set(CMAKE_CXX_COMPILER ${HOMEBREW_CLANGXX})
  else()
    message(WARNING "Homebrew Clang and Clang++ compilers not found. Using default compilers. 
                      OpenMP will be disabled.")
  endif()
endif()

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Set build flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -DLLVM_ENABLE_ASSERTIONS=ON")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# Find Python
set(Python_FIND_VIRTUALENV ONLY)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
message("Python3_EXECUTABLE " ${Python3_EXECUTABLE})

# # Find pybind11
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import site; print([s for s in site.getsitepackages() if 'site-packages' in s][0])" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
message("PYTHON_SITE_PACKAGES " ${PYTHON_SITE_PACKAGES})
set(pybind11_DIR ${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11)

message("pybind11_DIR " ${pybind11_DIR})
find_package(pybind11 REQUIRED)
message("pybind11_INCLUDE_DIRS " ${pybind11_INCLUDE_DIRS})

# *** BOOST SUPPORT: Find and include Boost ***
# find_package(Boost REQUIRED)
# message("Boost include dirs: ${Boost_INCLUDE_DIRS}")
# include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Try to find Boost (without REQUIRED, so we can handle its absence)
find_package(Boost)

if(Boost_FOUND)
  message("Boost include dirs: ${Boost_INCLUDE_DIRS}")
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
else()
  message(STATUS "Boost not found. Defining AMGCL_NO_BOOST to disable Boost-dependent features in AMGCL.")
  add_compile_definitions(AMGCL_NO_BOOST)
endif()

# Set include directories
set(DTCC_CPP_DIR ${CMAKE_SOURCE_DIR}/dtcc_core/cpp)
include_directories(${pybind11_INCLUDE_DIRS})
include_directories(${DTCC_CPP_DIR}/include)
include_directories(${DTCC_CPP_DIR}/include/citymodel)
include_directories(SYSTEM ${DTCC_CPP_DIR}/external)
include_directories(SYSTEM ${DTCC_CPP_DIR}/external/amgcl)

# Set installation directory
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

add_subdirectory(dtcc_core)

