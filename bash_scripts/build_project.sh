#!/bin/bash

# For an entire project

# MCMC Simulation lib
parent_dir="$(dirname -- "$(readlink -f -- "build_project.sh")")"
path_to_lattice_simulation_lib="$(dirname "$parent_dir")"

path_to_base_lib=${path_to_lattice_simulation_lib}

# Submodules
path_to_mcmc_simulation_lib="${path_to_lattice_simulation_lib}/external_submodules/MCMCSimulationLib/"
path_to_param_helper="${path_to_mcmc_simulation_lib}/external_submodules/ParamHelper/"

# Determine project_name, project_path and project_name
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_project_builder.sh"

# Build the project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_main_builder.sh"

# Compile the sample project
source "${path_to_mcmc_simulation_lib}/bash_scripts/generic_compiling.sh"

# Add bash_script for possibility to generate simulations
mkdir -p "$project_path/bash_scripts/"
source "${path_to_lattice_simulation_lib}/bash_scripts/write_project_related_build_simulation_sh_file.sh"

# Copy python scripts
cp -r "$path_to_lattice_simulation_lib/bash_scripts/python_scripts/" "${project_path}/"

# Copy jupyter notebooks
cp -r "$path_to_lattice_simulation_lib/bash_scripts/jupyter_notebooks/" "${project_path}/"

# Copy raw_transformer.py
mkdir -p "${project_path}/data/ONModelMetropolis/raw/"
cp "$path_to_lattice_simulation_lib/bash_scripts/raw_transformer/raw_transformer.py" "${project_path}/data/ONModelMetropolis/raw/raw_transformer.py"