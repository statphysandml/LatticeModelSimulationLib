pybind11_add_module(latticemodelsimulationlib_python LatticeModelSimulationLib_python.cpp)

# General compilation options- Already set
# set(RUN_WITH_PYTHON_BACKEND ON CACHE BOOL "Use Python backend to compute expectation values, etc.")

# Project specific compliation options (used by config.h.in)
set(CLUSTER_MODE "local" CACHE STRING "Mode to run the cluster.")
set(PYTHON_SCRIPTS_PATH "./../../python_scripts" CACHE STRING "Path to the potential python scripts evaluated by the python backend")

configure_file(config.h.in ./../../python/config.h @ONLY)

# Added so that 
find_library(ParamHelper NAMES libparamhelper.a PATHS ../build/ext/MCMCSimulationLib/ext/ParamHelper/src)
message("ParamHelper = ${ParamHelper}")
include_directories(../ext/MCMCSimulationLib/ext/ParamHelper/include/)

find_library(MCMCSimulationLib NAMES libmcmcsimulationlib.a PATHS ../build/ext/MCMCSimulationLib/src)
message("MCMCSimulationLib = ${MCMCSimulationLib}")
include_directories(../ext/MCMCSimulationLib/include/)

target_link_libraries(latticemodelsimulationlib_python PUBLIC latticemodelsimulationlib ${MCMCSimulationLib} ${ParamHelper})
set_target_properties(latticemodelsimulationlib_python PROPERTIES OUTPUT_NAME latticemodelimpl)
