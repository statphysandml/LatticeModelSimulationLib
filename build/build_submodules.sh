git submodule update --init --recursive

# Build MCMCSimulationLib
cd ../external_submodules/MCMCSimulationLib/build
source build.sh

# Navigate back to build directory
cd ../../../build/