add_library(latticemodelsimulationlib STATIC
  representations/links/u1.cpp
)

if (RUN_WITH_PYTHON_BACKEND)
  target_compile_definitions(latticemodelsimulationlib PUBLIC RUN_WITH_PYTHON_BACKEND)
  target_link_libraries(latticemodelsimulationlib PRIVATE ${target_link_libraries_python_appendix})
endif()

target_include_directories(latticemodelsimulationlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
