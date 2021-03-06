#!/bin/bash

# For small simulations

# Lattice Model Simulation Lib
if [ -z ${path_to_lattice_simulation_lib=+x} ]; then
  echo "Building a simulation without a project is currently not working because of a potentially wrong python_scripts path. But simulations can be generated with the bash script in the project_path/bash_scripts/ directory."
  exit
  parent_dir="$(dirname -- "$(readlink -f -- "build_simulation.sh")")"
  path_to_lattice_simulation_lib="$(dirname "$parent_dir")"
fi

path_to_base_lib=${path_to_lattice_simulation_lib}

# Submodules
path_to_mcmc_simulation_lib="${path_to_lattice_simulation_lib}/external_submodules/MCMCSimulationLib/"
path_to_param_helper="${path_to_mcmc_simulation_lib}/external_submodules/ParamHelper/"

# Determine project_name, project_path and project_type
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_simulation_builder.sh"

# Build the project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_main_builder.sh"

# Compile the sample project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_compiling.sh"
