#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <iostream>

#include "config.h"

// #include "LatticeModelSimulationLib/LatticeModelSimulationLib.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

// namespace LatticeModelSimulationLib {

#include <lattice_model_impl/representations/link_header.hpp>

#include <mcmc_simulation/header.hpp>
#include <execution/executer.hpp>

#include <lattice_model_impl/update_dynamics/update_dynamics_header.hpp>
#include <lattice_model_impl/mcmc_update/mcmc_update_header.hpp>

#include <lattice_model_impl/site/site_header.hpp>
#include <lattice_model_impl/lattice/lattice_header.hpp>
#include <lattice_model_impl/link_lattice/link_lattice_header.hpp>

using namespace lm_impl;




void init_lattice_update_formalisms(py::module &m)
{
    py::class_<update_dynamics::SequentialUpdateParameters>(m, "SequentialUpdateParameters")
        .def(py::init());    
}


void init_lattice_models(py::module &m)
{
    py::class_<lattice_system::IsingModelParameters>(m, "IsingModelParameters")
        .def(py::init<double, double, double>(), "beta"_a, "J"_a, "h"_a)
        .def_readonly("beta", &lattice_system::IsingModelParameters::beta)
        .def_readonly("J", &lattice_system::IsingModelParameters::J)
        .def_readonly("h", &lattice_system::IsingModelParameters::h);
}


void init_ising_model(py::module &m)
{
    // MCMC Update Formalism
    typedef mcmc_update::MetropolisUpdateParameters<lattice_system::IsingModelParameters, lattice_system::IsingModelSampler> IsingModelMetropolisParameters;
    py::class_<IsingModelMetropolisParameters>(m, "IsingModelMetropolisParameters")
        .def(py::init<double>(), "eps"_a=0);
        // .def_readonly("eps", &IsingModelMetropolisParameters::eps);
    
    // Lattice System Parameters
    typedef lattice_system::LatticeParameters<double, lattice_system::IsingModelParameters, IsingModelMetropolisParameters, update_dynamics::SequentialUpdateParameters> IsingModelMetropolisLatticeParameters;
    py::class_<IsingModelMetropolisLatticeParameters>(m, "IsingModelMetropolisLatticeParameters")
        .def_static("generate_parameters",
            py::overload_cast<
                lattice_system::IsingModelParameters&,
                IsingModelMetropolisParameters&,
                update_dynamics::SequentialUpdateParameters&,
                const std::vector<int>&,
                const std::string&,
                const std::vector<std::string>&>(&IsingModelMetropolisLatticeParameters::generate_parameters),
            "model_parameters"_a, "mcmc_update_parameters"_a, "latice_update_parameters"_a, "dimensions"_a=std::vector<int>{4, 4}, "lattice_action_type"_a="nearest_neighbour", "measures"_a=std::vector<std::string>{})
        .def_readonly("dimensions", &IsingModelMetropolisLatticeParameters::dimensions)
        .def_readonly("lattice_action_type", &IsingModelMetropolisLatticeParameters::lattice_action_type);
        // .def("update_measures", &IsingModelMetropolisLatticeParameters::update_measures, "measures"_a);

    // Lattice System
    typedef lattice_system::LatticeSystem<double, lattice_system::IsingModelParameters, IsingModelMetropolisParameters, update_dynamics::SequentialUpdateParameters> IsingModelMetropolisLattice;
    py::class_<IsingModelMetropolisLattice>(m, "IsingModelMetropolisLattice")
        .def(py::init<IsingModelMetropolisLatticeParameters&>(), "IsingModelMetropolisLatticeParameters"_a)
        .def_readonly("lattice", &IsingModelMetropolisLattice::lattice)
        .def("update_step", &IsingModelMetropolisLattice::update_step, "measure_interval"_a)
        .def("initialize", &IsingModelMetropolisLattice::initialize, "starting_mode"_a)
        .def("measure", &IsingModelMetropolisLattice::measure)
        .def("measure_names", &IsingModelMetropolisLattice::measure_names);

    // Simulation Parameters

    // -> Equilibrium Time
    typedef mcmc::simulation::SimulationParameters<IsingModelMetropolisLatticeParameters, mcmc::execution::EquilibriumTimeParameters> IsingModelMetropolisEquilibriumTimeParameters;
    py::class_<IsingModelMetropolisEquilibriumTimeParameters>(m, "IsingModelMetropolisEquilibriumTimeParameters")
        .def_static("generate_simulation",
            py::overload_cast<
                IsingModelMetropolisLatticeParameters&,
                mcmc::execution::EquilibriumTimeParameters&,
                const std::string&,
                const std::string&,
                const std::string&,
                const double&,
                const double&,
                const int&>(&IsingModelMetropolisEquilibriumTimeParameters::generate_simulation),
            "systembase_parameters"_a, "execution_parameters"_a, "rel_data_path"_a,
            "running_parameter_kind"_a="None", "running_parameter"_a="None",
            "rp_minimum"_a=0.0, "rp_maximum"_a=0.0, "rp_number"_a=0.0)
        .def_readonly("rp_minimum", &IsingModelMetropolisEquilibriumTimeParameters::rp_minimum)
        .def_readonly("rp_maximum", &IsingModelMetropolisEquilibriumTimeParameters::rp_maximum)
        .def_readonly("rp_number", &IsingModelMetropolisEquilibriumTimeParameters::rp_number);

    // -> Correlation Time
    typedef mcmc::simulation::SimulationParameters<IsingModelMetropolisLatticeParameters, mcmc::execution::CorrelationTimeParameters> IsingModelMetropolisCorrelationTimeParameters;
    py::class_<IsingModelMetropolisCorrelationTimeParameters>(m, "IsingModelMetropolisCorrelationTimeParameters")
        .def_static("generate_simulation",
            py::overload_cast<
                IsingModelMetropolisLatticeParameters&,
                mcmc::execution::CorrelationTimeParameters&,
                const std::string&,
                const std::string&,
                const std::string&,
                const double&,
                const double&,
                const int&>(&IsingModelMetropolisCorrelationTimeParameters::generate_simulation),
            "systembase_parameters"_a, "execution_parameters"_a, "rel_data_path"_a,
            "running_parameter_kind"_a="None", "running_parameter"_a="None",
            "rp_minimum"_a=0.0, "rp_maximum"_a=0.0, "rp_number"_a=0.0)
        .def_readonly("rp_minimum", &IsingModelMetropolisCorrelationTimeParameters::rp_minimum)
        .def_readonly("rp_maximum", &IsingModelMetropolisCorrelationTimeParameters::rp_maximum)
        .def_readonly("rp_number", &IsingModelMetropolisCorrelationTimeParameters::rp_number);

    // -> Expectation Value
    typedef mcmc::simulation::SimulationParameters<IsingModelMetropolisLatticeParameters, mcmc::execution::ExpectationValueParameters> IsingModelMetropolisExpectationValueParameters;
    py::class_<IsingModelMetropolisExpectationValueParameters>(m, "IsingModelMetropolisExpectationValueParameters")
        .def_static("generate_simulation",
            py::overload_cast<
                IsingModelMetropolisLatticeParameters&,
                mcmc::execution::ExpectationValueParameters&,
                const std::string&,
                const std::string&,
                const std::string&,
                const double&,
                const double&,
                const int&>(&IsingModelMetropolisExpectationValueParameters::generate_simulation),
            "systembase_parameters"_a, "execution_parameters"_a, "rel_data_path"_a,
            "running_parameter_kind"_a="None", "running_parameter"_a="None",
            "rp_minimum"_a=0.0, "rp_maximum"_a=0.0, "rp_number"_a=0.0)
        .def_readonly("rp_minimum", &IsingModelMetropolisExpectationValueParameters::rp_minimum)
        .def_readonly("rp_maximum", &IsingModelMetropolisExpectationValueParameters::rp_maximum)
        .def_readonly("rp_number", &IsingModelMetropolisExpectationValueParameters::rp_number);
    
    // Simulation

    // -> Equilibrium Time
    typedef mcmc::simulation::Simulation<IsingModelMetropolisLatticeParameters, mcmc::execution::EquilibriumTimeParameters> IsingModelMetropolisEquilibriumTime;
    py::class_<IsingModelMetropolisEquilibriumTime>(m, "IsingModelMetropolisEquilibriumTime")
        .def(py::init<IsingModelMetropolisEquilibriumTimeParameters&>())
        .def("run", &IsingModelMetropolisEquilibriumTime::run)
        .def("evaluate", &IsingModelMetropolisEquilibriumTime::evaluate, "rel_results_dir"_a, "sim_root_dir"_a);

    // -> Correlation Time
    typedef mcmc::simulation::Simulation<IsingModelMetropolisLatticeParameters, mcmc::execution::CorrelationTimeParameters> IsingModelMetropolisCorrelationTime;
    py::class_<IsingModelMetropolisCorrelationTime>(m, "IsingModelMetropolisCorrelationTime")
        .def(py::init<IsingModelMetropolisCorrelationTimeParameters&>())
        .def("run", &IsingModelMetropolisCorrelationTime::run)
        .def("evaluate", &IsingModelMetropolisCorrelationTime::evaluate, "rel_results_dir"_a, "sim_root_dir"_a);

    // -> Expectation Value
    typedef mcmc::simulation::Simulation<IsingModelMetropolisLatticeParameters, mcmc::execution::ExpectationValueParameters> IsingModelMetropolisExpectationValue;
    py::class_<IsingModelMetropolisExpectationValue>(m, "IsingModelMetropolisExpectationValue")
        .def(py::init<IsingModelMetropolisExpectationValueParameters&>())
        .def("run", &IsingModelMetropolisExpectationValue::run)
        .def("evaluate", &IsingModelMetropolisExpectationValue::evaluate, "rel_results_dir"_a, "sim_root_dir"_a);
}

PYBIND11_MODULE(latticemodelimpl, m)
{
    // init_functions(m);
    // init_execution_modes(m);
    init_lattice_update_formalisms(m);
    init_lattice_models(m);

    init_ising_model(m);

    m.doc() = "Python Bindings for LatticeModelImpl";
} // namespace LatticeModelSimulationLib
