LatticeModelSimulationLib
=========================

<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/statphysandml/LatticeModelSimulationLib/CI)](https://github.com/statphysandml/LatticeModelSimulationLib/actions?query=workflow%3ACI)
[![PyPI Release](https://img.shields.io/pypi/v/LatticeModelSimulationLib.svg)](https://pypi.org/project/LatticeModelSimulationLib)
[![Documentation Status](https://readthedocs.org/projects/LatticeModelSimulationLib/badge/)](https://LatticeModelSimulationLib.readthedocs.io/)
[![codecov](https://codecov.io/gh/statphysandml/LatticeModelSimulationLib/branch/main/graph/badge.svg)](https://codecov.io/gh/statphysandml/LatticeModelSimulationLib)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=statphysandml_LatticeModelSimulationLib&metric=alert_status)](https://sonarcloud.io/dashboard?id=statphysandml_LatticeModelSimulationLib) -->

LatticeModelSimulationLib is C++ library which is built upon the MCMCSimulationLib. It implements different well-known Monte Carlo algorithms and interesting models in statistical physics. In the long run, the library should be a collection for all common used models. Switching between different models and different algorithms is due to the modular structure very easy.

The evaluation of the simulation data takes place in Python. The Python modules allow a convenient computation of observables and further interesting quantities. In addition, the Python modules can be used to transform the data into a PyTorch dataset and to train a machine learning algorithm on the data.

The MCMCSimulationLib is integrated into the project as a submodule in the ext/ directory.

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

Most of the implemented models are presented in more detail together with examples in the example/ directory.

More models and algorithms will follow. Feel also free to contribute to this library!

In general, the library provides code for the implementation of site and lattice models and can handle arbitrary dimensions.

It is straightforward to implement a new model that makes use of the existing algorithms and implementation. A possible approach is to copy the code of an exsiting model and adapt it to your new model.

Furthermore, a python wrapper will be implemented to provide an easier access for using the library.

Prerequisites
-------------

Building LatticeModelSimulationLib requires the following software installed:

* A C++14-compliant compiler
* CMake `>= 3.15`
* Cookiecutter, e.g. by doing pip install cookiecutter, for auto-generating new projects
* Doxygen (optional, documentation building is skipped if missing) (not present, yet)
* Python `>= 3.6` for building Python bindings and for running the evaluation scripts of the library (will follow)
* Boost `>= 1.67` (https://www.boost.org/)

as well as the following python dependencies for a successful usage of the evaluation of numerical results:

* [MCMCEvaluationLib](https://github.com/statphysandml/MCMCEvaluationLib)
* [pystatsplottools](https://github.com/statphysandml/pystatplottools)


Building LatticeModelSimulationLib
----------------------------------

The following sequence of commands builds LatticeModelSimulationLib. If you use a virtual envirnonment, it is important that it is activated for the building process to find the correct python version. The sequence assumes that your current working directory is the top-level directory
of the freshly cloned repository:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

The build process can be customized with the following CMake variables,
which can be set by adding `-D<var>={ON, OFF} / {other}` to the `cmake` call:

* `RUN_WITH_PYTHON_BACKEND`: Enable the computation of numerical results in Python (default: `ON`). Some features, like the computation of the correlation time do not work if this is disabled.
* `BUILD_TESTING`: Enable building of the test suite (default: `ON`)
* `BUILD_DOCS`: Enable building the documentation (default: `ON`)
<!-- * `BUILD_PYTHON_BINDINGS`: Enable building the Python bindings (default: `ON`) -->

The following additional CMake variables are only required if the library is used to submit jobs on a cluster:

* `VIRTUAL_ENV`: Name of the virtual environment to be acivated before execution of the script
* `CONDA_ACTIVATE_PATH`: Path conda activate (for example: "~/miniconda3/bin/activate")

Note that the library is built in almost the exact same way as the MCMCSimulationLib (https://github.com/statphysandml/MCMCSimulationLib). For the MCMCSimulationLib there exists a more thorough getting started guide (https://github.com/statphysandml/MCMCSimulationLib/blob/master/doc/getting_started.md) that also explains how the virtual environment for Python can be set up. It also contains a short introduction to the LatticeModelSimulationLib. But it is important to note that building the MCMCSimulationLib is not necesseary for this project since the library is integrated here as a submodule.

<!-- # Documentation

LatticeModelSimulationLib provides a Sphinx-based documentation, that can
be browsed [online at readthedocs.org](https://LatticeModelSimulationLib.readthedocs.io). -->

Examples
--------

Several examples can be found in a separate repository: https://github.com/statphysandml/LatticeModelImplementations. The respository is also included as a submodule in the examples/ directory.

After a navigation into the top-level directory, the following sequence of commands builds the respective example. Note again that if you use a virtual envirnonment, it is important that it is activate for the building process to find the correct python version.

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

The build process can be again customized with the following CMake variables,
which can be set by adding `-D<var>={ON, OFF} / {other}` to the `cmake` call:

* `CLUSTER_MODE`: Indicates whether the code is executed on the cluster (`on_cluster`) or locally (`local`) (default: `local`). The "local" option can be used to test whether the library prepares a computation on a cluster correctly. In this case, the simulation runs locally on your machine. The option can be changed to "on_cluster". In this case the jobs are sent to the cluster. There are two functions "prepare_execution_on_cpu_cluster" and "run_execution_on_cpu_cluster" that take care of this. The functions can be found in the file ext/MCMCSimulationLib/src/execution/execution.cpp and need to be adapted according to the used cluster. More details on how to execute a simulation on the cluster can be found in the main.cpp file of the SimulateAndExecute example (https://github.com/statphysandml/MCMCSimulationLib/blob/master/examples/SimulateAndExecute/src/main.cpp) or in the main.cpp file of the template project (see Template Project).
* `PYTHON_SCRIPTS_PATH`: Path to a directory including additional python files for a possible execution of code of custom written functions and modules. (default: `./python_scripts`). The path is appended by the programming to sys.path. It needs to be defined relative to the project root path.

Template Project
----------------

A new project demonstrating most of the functionalities of the library based on a simulation of the O(n) model with a Metropolis algorithm can be generated by running:

```python
python generate_application.py [-h] [-o OUTPUT_DIR] [--config-file CONFIG_FILE] [-n PROJECT_NAME]
```

where the arguments have the following meanings:

* -o, --output-dir: Where to output the generated project dir into.
* --config-file, User configuration file for overwriting default parameters of cookiecutter.json file in the application_wrapper/ directory (default: None)
* -n, --project_name: The name of the project (default: "my-mcmc-simulation-project")

Instruction for building the project can be found in a respective generated README.md file or in the previous examples section (which works equivalently).

The generated src/main.cpp file contains detailed instructions on different ways to execute the simulation. The project can be used as a starting point for your own simulations.

The template project demonstrates also some ways to work with the resulting data in Python. In particular a jupyter notebook is generated in the jupyter_notebooks/ directory that can be executed after a successfull simulation of the O(n) model. This also includes a generation of a .pt (PyTorch data file). This makes the training of a machine learning algorithms with the sampled configurations pretty easy.

Usage
-----

cmake (CMakeLists.txt) - The CMakeLists.txt files of the template project and the examples/ directory contain examples for a possible integration of the library with CMake.

Further comments
----------------

To keep track of changes in the submodules when pulling updates of the respository, we strongly recommend for an update on your local machine to use

```bash
git pull --recurse-submodules
```

otherwise the code might not run properly after an update of the library.


Projects using the MCMCSimulationLib library
-------------------------------------------------------

- LatticeModelImplementations (https://github.com/statphysandml/LatticeModelImplementations)

Support and Development
----------------------

The project was generated with the help of the [cookiecutter-cpp-project](https://github.com/ssciwr/cookiecutter-cpp-project) of the [Scientific Software Center, IWR, Heidelberg University](https://ssc.iwr.uni-heidelberg.de/).

For bug reports/suggestions/complaints please file an issue on GitHub.

Or start a discussion on our mailing list: statphysandml@thphys.uni-heidelberg.de
