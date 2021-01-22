LatticeModelSimulationLib
=================

LatticeModelSimulationLib is C++ library which is build upon the MCMCSimulationLib. It implements different well-known Monte Carlo algorithms and interesting models in statistical physics. In the long run, the library should be a collection for all common used models. Switching between different models and different algorithms is due to the modular structure very easy.

The evaluation of the simulation data takes place in Python. The Python modules allow a convenient computation of observables and further interesting quantities. In addition, the Python modules can be used to transform the data into a PyTorch dataset and to train a machine learning algorithm on the data.

The MCMCSimulationLib is integrated into the project as a submodule in the external_submodule/ directory.

The library implements the following Monte Carlo algorithms:

- the Metropolis algorithm
- the Hybrid/Hamiltonian Monte Carlo algorithm
- Langevin dynamics
- Complex Langevin dynamics

the following lattice models:

- the Ising model
- the classical XY model
- the O(n) model
- the U(1) model
- the SU(2) model
- the complex ON model
- the complex XY model

and the following site models:

- the complex Gaussian model
- the complex cubic model
- the complex and real polynomial/quartic model

Most of the implemented models are presented in more detail together with examples here: https://github.com/statphysandml/LatticeModelImplementations.

More models and algorithms will follow. Feel also free to contribute to this library!

In general, the library provides code for the implementation of site and lattice models and can handle arbitrary dimensions.

It is straightforward to implement a new model that makes use of the existing algorithms and implementation. A possible approach is to copy the code of an exsiting model and adapt it to your new model.

Build
-----

The library is build in almost the exact same way as the MCMCSimulationLib (https://github.com/statphysandml/MCMCSimulationLib). Therefore, the most of instruction is taken from there:

Certain configuration parameters need to be defined for an execution of the C++ program. They need to be defined in a config.sh file in the build/ directory and in a project_config.sh file in the bash_scripts/ directory. The config.sh file contains all important information for the integration of Python into the program. So far, the library can only be used with a virtual environment. The config_template.sh file is a template where all necessary parameters are defined:
```bash
path_to_python3="~/.miniconda3/envs/virtual_env/" # (optional)
virtual_env="virtual_env" # (optional)
python_version="3.7" # (optional)
path_to_conda_activate="~/.miniconda3/bin/activate" # (optional)

# Optional

# path_to_boost="/opt/boost_1_70_0" # can be defined if CMake cannot find boost

# Only necessary, if GPU is used - adapt this to your architecture
# path_to_cuda="/opt/cuda-10.1" # can be defined if CMake cannot find cuda
# nvcc_flag_gencode_arch=compute_75
# nvcc_flag_gencode_code=sm_75
```
The parameters need to be adapted to the virtual environment of the underlying system. The parameter path_to_boost can be used if CMake cannot find boost. The same holds for CUDA. However, no functionalities of CUDA are used within the library at the moment, so nothing needs to be done here, for now.

Using the library without any Python components is possible by leaving the config.sh file empty.

Having a config.sh file in the build directory, the library can be build with
```bash
cd build
bash build.sh
```
The file project_config.sh defines additional project-dependent parameters. They are used when a new project is generated with the build_project.sh file. There also exists a template file (project_config_template.sh) which contains all necessary parameters:
```bash
cluster_mode="local"

# (optional - default=local) local/on_cluster - Can be adapted temporarily by
# adding -DCLUSTER_MODE=".." to the cmake command
# - "local" = for testing - does not actually start the code on the cluster but locally
# and performs all the necessary preparation
# - "on_cluster" = for the actual execution on a cluster

# python_modules_path="./python_scripts"

# (optional - default="./python_scripts" for projects and "./../python_scripts/" for
# simulations.) for a possible execution code of custom written functions and modules.
# The directory "python_modules_path" is added to sys.path by the program.
# The path needs to be defined relative to the project root path)
```

The first parameter "cluster_mode" indicates whether the algorithms are started on the cluster (on_cluster) or locally (local).

#### cluster_mode = local

The "local" option can be used to test whether the library prepares a computation on a cluster correctly. In this case, the simulation runs locally on your machine.

#### cluster_mode = on_cluster

The option can be changed to "on_cluster". In this case the jobs are sent to the cluster. There are two functions "prepare_execution_on_cpu_cluster" and "run_execution_on_cpu_cluster" that take care of this. The functions can be found in the file src/execution/execution.cpp of the MCMCSimulationLib (https://github.com/statphysandml/MCMCSimulationLib/blob/master/src/execution/executer.cpp) and need to be adapted according to the used cluster.

Both options can be updated by building the CMake files of the executable again. You might change the cluster mode, for example, by executing
```bash
cmake ../cmake/ -DCMAKE_BUILD_TYPE=Release -DCLUSTER_MODE=on_cluster
make -j4
```

in the release directory of your project. More details on how to execute a simulation on the cluster can be found in the main.cpp file of the SimulateAndExecute example (https://github.com/statphysandml/MCMCSimulationLib/blob/master/examples/simulations/SimulateAndExecute/main.cpp) or in the main.cpp file of a template project (see Template Project).

Examples
--------

Several examples can be found in a separate project: https://github.com/statphysandml/LatticeModelImplementations.

Template Project
----------------

A new project that demonstrates most of the functionalities of the library based on a the simulation of the O(n) model with a Metropolis algorithm can be generated in an arbitrary directory with

```bash
cd bash_scripts
bash build_project.sh
```

A project can be used as a template for your own simulation. The main.cpp contains detailed instructions on different ways to execute the simulation and provides a good starting point. The project itself has a directory bash_scripts/ with an additional build_simulation.sh bash script. This can be used to generate a template simulation for a further executable which uses code of your project. These smaller simulation projects are very useful since they enable the generation of exeuctables for several different simulations of your projects. This enables a rerunning of a simulation at a later point and simplifies running code on a cluster. For example, one might use the following workflow:

- Write all important code in the include/ and the src/ directory of your project.
- Generate a new simulation with the build_simulation.sh file.
- Adapt the main.cpp file in your simulation for an execution of a MCMC simulation.
- Transfer the code to your cluster and rerun build_simulation.sh. Enter the name of your simulation. The bash script recognizes that the simulation already exists and just adapts the CMakeLists.txt file according to the settings in your project_config.sh file.
- Execute the code with the mcmc::exeuction::execute<>() function to submit the simulation to the cluster.

Both, the build_project.sh and the build_simulation.sh bash script recognize if a project or a simulation has been built before. If these scripts are executed with the same input, only the CMakeLists.txt file is adapted according to the settings in the project_config.sh file and with respect to the executables of the local libraries (which can be found in the external_submodules/ directory). Note that it is sometimes necessary to delete the files in the debug or release that are generated from CMake if parameters, for example, in the config files, have been changed.

The template project demonstrates also some ways to work with the resulting data in Python. In particular a jupyter notebook is generated in the jupyter_notebooks/ directory that can be executed after a successfull simulation of the O(n) model. This also includes a generation of a .pt (PyTorch data file). This makes the training of a machine learning algorithms with the sampled configurations pretty easy.

Usage
-----

cmake (CMakeLists.txt) - The CMakeLists.txt files of the template project contains an example for a possible integration of the library with CMake.

Dependencies
------------

- boost (https://www.boost.org/)

Projects using the MCMCSimulationLib library
-------------------------------------------------------

- LatticeModelImplementations (https://github.com/statphysandml/LatticeModelImplementations)

Support and development
----------------------

For bug reports/suggestions/complaints please file an issue on GitHub.

Or start a discussion on our mailing list: statphysandml@thphys.uni-heidelberg.de
