set(BINDINGS _dtcc_builder)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)

pybind11_add_module(${BINDINGS} ${SRC_DIR}/dtcc_builder.cpp)
target_compile_features(${BINDINGS} PRIVATE cxx_std_17)
set_target_properties(${BINDINGS} PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(${BINDINGS} PRIVATE ${SRC_DIR}/include)
target_include_directories(${BINDINGS} SYSTEM PRIVATE ${SRC_DIR}/external)

if(DTCC_ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(${BINDINGS} PRIVATE /W4 /permissive-)
  else()
    target_compile_options(${BINDINGS} PRIVATE
      -Wall -Wextra -Wpedantic -Wno-reorder -Wno-unused-parameter)
  endif()
endif()

find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${BINDINGS} PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not available; continuing without it")
endif()

if(DTCC_USE_TRIANGLE)
  # The bundled triangle.h contains the implementation, not just declarations.
  target_include_directories(${BINDINGS} SYSTEM PRIVATE ${TRIANGLE_INCLUDE_DIR})
  target_compile_definitions(${BINDINGS} PRIVATE DTCC_HAVE_TRIANGLE TRILIBRARY ANSI_DECLARATORS)
  if(WIN32)
    target_compile_definitions(${BINDINGS} PRIVATE NO_TIMER)
  endif()
endif()

if(APPLE)
  target_link_options(${BINDINGS} PRIVATE "-Wl,-undefined,error")
endif()

install(TARGETS ${BINDINGS}
  LIBRARY DESTINATION dtcc_core/builder
  ARCHIVE DESTINATION dtcc_core/builder
  RUNTIME DESTINATION dtcc_core/builder
)
