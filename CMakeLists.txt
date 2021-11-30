cmake_minimum_required(VERSION 3.15)

# Set a name and a version number for your project:
project(LatticeModelSimulationLib VERSION 0.0.1 LANGUAGES CXX)

# Initialize some default paths
include(GNUInstallDirs)

# Define the minimum C++ standard that is required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Compilation options
set(BUILD_PYTHON_BINDINGS ON CACHE BOOL "Enable building of Python bindings")
set(RUN_WITH_PYTHON_BACKEND ON CACHE BOOL "Use Python backend to compute expectation values, etc.")
set(BUILD_DOCS ON CACHE BOOL "Enable building of documentation")

# Additional compliation options to be able to activate the environment before the actual simulation in a generated bash script.
# These are only required for runs with a queue-script. (used by include/execution/config.h.in)
set(VIRTUAL_ENV "None" CACHE STRING "Name of the virtual environment")
set(CONDA_ACTIVATE_PATH "None" CACHE STRING "Path to conda activate path-to-conda-activate/bin/activate")

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/ext/MCMCSimulationLib/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()


configure_file(./ext/MCMCSimulationLib/include/execution/config.h.in ../ext/MCMCSimulationLib/include/execution/config.h @ONLY)
# configure_file(./include/execution/config.h.in ../include/execution/config.h @ONLY)

# Boost
FIND_PACKAGE( Boost 1.67 REQUIRED COMPONENTS filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("Boost = ${Boost_INCLUDE_DIRS}")
endif()

# Compilation options

# compile the external submodule
add_subdirectory(ext/MCMCSimulationLib)
include_directories(ext/MCMCSimulationLib/include)

# compile the library
add_subdirectory(src)


# compile the tests
include(CTest)
if(BUILD_TESTING)
  # add_subdirectory(ext/Catch2)
  include(./ext/MCMCSimulationLib/ext/ParamHelper/ext/Catch2/contrib/Catch.cmake)
  add_subdirectory(tests)
endif()

if(BUILD_DOCS)
  # Add the documentation
  add_subdirectory(doc)
endif()


# Add an alias target for use if this project is included as a subproject in another project
add_library(latticemodelsimulationlib::latticemodelsimulationlib ALIAS latticemodelsimulationlib)

# Install targets and configuration
install(
  TARGETS latticemodelsimulationlib
  EXPORT latticemodelsimulationlib-config
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
  EXPORT latticemodelsimulationlib-config
  NAMESPACE latticemodelsimulationlib::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/latticemodelsimulationlib
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# This prints a summary of found dependencies
include(FeatureSummary)
feature_summary(WHAT ALL)
