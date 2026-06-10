cmake_minimum_required(VERSION 3.18)

set(BINDINGS _dtcc_builder)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)

pybind11_add_module(${BINDINGS} ${SRC_DIR}/dtcc_builder.cpp)

target_include_directories(${BINDINGS}
  PRIVATE
    ${DTCC_CORE_PUBLIC_INCLUDE_DIRS}
)

if(DTCC_CORE_SYSTEM_INCLUDE_DIRS)
  target_include_directories(${BINDINGS} SYSTEM PRIVATE ${DTCC_CORE_SYSTEM_INCLUDE_DIRS})
endif()

if(DTCC_ENABLE_WARNINGS AND NOT MSVC)
  target_compile_options(${BINDINGS} PRIVATE -Wno-unused-parameter)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${BINDINGS} PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not available; continuing without it")
endif()

if(DTCC_CORE_LINK_LIBS)
  target_link_libraries(${BINDINGS} PRIVATE ${DTCC_CORE_LINK_LIBS})
endif()

if(DTCC_HAVE_BOOST)
  target_include_directories(${BINDINGS} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
else()
  target_compile_definitions(${BINDINGS} PRIVATE AMGCL_NO_BOOST)
endif()

if(DTCC_HAVE_TRIANGLE)
  target_include_directories(${BINDINGS} SYSTEM PRIVATE ${TRIANGLE_INCLUDE_DIR})
  target_link_libraries(${BINDINGS} PRIVATE ${TRIANGLE_LIBRARY})
  target_compile_definitions(${BINDINGS} PRIVATE DTCC_HAVE_TRIANGLE TRILIBRARY ANSI_DECLARATORS)
endif()

if(WIN32)
  target_compile_definitions(${BINDINGS} PRIVATE NO_TIMER)
endif()

if(APPLE)
  target_link_libraries(${BINDINGS} PRIVATE c++)
  target_link_options(${BINDINGS} PRIVATE "-Wl,-undefined,error")
endif()

install(TARGETS ${BINDINGS}
  LIBRARY DESTINATION dtcc_core/builder
  ARCHIVE DESTINATION dtcc_core/builder
  RUNTIME DESTINATION dtcc_core/builder
)

if(NOT SKBUILD)
  install(TARGETS ${BINDINGS}
    LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/dtcc_core/builder
    ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/dtcc_core/builder
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/dtcc_core/builder
  )
endif()
