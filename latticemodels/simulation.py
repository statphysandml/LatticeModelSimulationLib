from pystatplottools.pdf_env.loading_figure_mode import loading_figure_mode
fma, plt = loading_figure_mode(develop=True)

import latticemodelimpl
import numpy as np
import pandas as pd

from latticemodels.ising_model import IsingModel
from latticemodels.mcmc.evaluation_module import EvaluationModule


""" -> split into two? On for custom simulations and one for simulations with the C++ functionality as a Wrapper.
 Allows the definition of a further BaseClass for lattice models.
 -> IsingModelSimulation, inheriting from SimulationModelBase -> for non custom simulations
 -> IsingModel, inheriting from EvaluationModule -> for custom simulations """

""" Allows for custom simulations -> no access to mode simulation tools! only to the one of mcmc simulation lib by providing directly the necessary data """


class LatticeSimulation(EvaluationModule):
    def __init__(self, model, default_measures=[], sim_base_dir=None,
                 rel_data_path=None,  # -> sim_base_dir + "/" + rel_data_path
                 rel_results_path=None,  # -> sim_base_dir + "/" + rel_results_path
                 running_parameter_kind=None,
                 running_parameter=None,
                 rp_keys=None, # -> either rp_keys or rp_minimum, etc.
                 rp_minimum=0.0,
                 rp_maximum=0.0,
                 rp_number=0):
        super().__init__(sim_base_dir=sim_base_dir, rel_data_path=rel_data_path, rel_results_path=rel_results_path,
                         running_parameter_kind=running_parameter_kind, running_parameter=running_parameter,
                         rp_keys=rp_keys, rp_minimum=rp_minimum, rp_maximum=rp_maximum, rp_number=rp_number) 

        self.model = model
        self.mcmc_system = None

        self.measurements = {}

    def initialize_lattice(self, rp_val=None, starting_mode="hot"):
        if rp_val is not None:
            setattr(self.model, self.running_parameter, rp_val)

        self.model.initialize_parameters()

        self.mcmc_system = self.model.generate_mcmc_system()
        self.mcmc_system.initialize(starting_mode=starting_mode)
    
    def update(self, n_steps):
        self.mcmc_system.update_step(n_steps)

    @property
    def mcmc_model(self):
        return self.mcmc_system.lattice

    def get_lattice(self):
        pass

    def at(self):
        pass

    def measure(self, running_parameter="default"):
        measurements = self.mcmc_system.measure()
        rp_val = getattr(self.model, self.running_parameter)
        assert self.mcmc_system.measure_names() == self.model.measures,\
            "Initialized measures and actually measured measurements do not coincide."
        for measurement, measure_name in zip(measurements, self.model.measures):
            self.measurements[rp_val][measure_name].append(measurement)

    def initialize_measurements(self, measures):
        self.model.set_measures(measures=measures)
        self.measurements = {rp_val: {measure: [] for measure in self.model.measures} for rp_val in self.rp_keys}

    def measurements_to_dataframe(self, complex_number_format="complex", transformer=None, transform=False, transformer_path=None):
        from mcmctools.loading.loading import ConfigurationLoader
        n_measurements = len(self.measurements[self.rp_keys[0]][self.model.measures[0]])
        data = ConfigurationLoader.process_mcmc_configurations(
            data=[pd.DataFrame({**item, self.running_parameter.capitalize(): [key] * n_measurements}) for key, item in self.measurements.items()],
            running_parameter=self.running_parameter, complex_number_format=complex_number_format, transformer=transformer, transform=transform, transformer_path=transformer_path)
        return data

    def measurements_to_file(self):
        pass

    def run_equilibrium_time_simulation(self, measure, sample_size, number_of_steps):
        self.initialize_measurements(measures=[measure])

        # Possibility to define a custom __iter__ class - or several for each mode...
        for rp_val in self.rp_keys:
            starting_mode = "hot"
            for m in range(2 * sample_size):
                self.initialize_lattice(starting_mode=starting_mode, rp_val=rp_val)
                self.measure()
                for n in range(number_of_steps - 1):
                    self.update(n_steps=1)
                    self.measure()

                if starting_mode == "hot":
                    starting_mode = "cold"
                else:
                    starting_mode = "hot"

    def run_correlation_time_simulation(self, measure, minimum_sample_size, maximum_correlation_time, start_measuring):
        self.initialize_measurements(measures=[measure])

        # Possibility to define a custom __iter__ class - or several for each mode...
        for rp_val in self.rp_keys:
            for m in range(minimum_sample_size):
                self.initialize_lattice(starting_mode="hot", rp_val=rp_val)
                self.update(n_steps=start_measuring)
                self.measure()
                for n in range(maximum_correlation_time - 1):
                    self.update(n_steps=1)
                    self.measure()

    def run_expectation_value_simulation(self, measures, n_measurements, n_steps_equilibrium, n_steps_autocorrelation, starting_mode="hot"):
        self.initialize_measurements(measures=measures)

        # Possibility to define a custom __iter__ class - or several for each mode...
        for rp_val in self.rp_keys:
            self.initialize_lattice(starting_mode=starting_mode, rp_val=rp_val)
            self.update(n_steps=n_steps_equilibrium)
            self.measure()
            for n in range(n_measurements - 1):
                self.update(n_steps=n_steps_autocorrelation)
                self.measure()


if __name__ == "__main__":
    # test_LatticeModelSimulationLib()
    # Overwrite data directory?
    ising_model = IsingModel(beta=0.1, J=1.0, h=0.0, dimensions=[4, 4])

    simulation = LatticeSimulation(model=ising_model,
                             rel_data_path="./data/Test/",
                             rel_results_path="./data/Test/results/",
                             # running_parameter_kind="model_params",
                             running_parameter="beta",
                             rp_keys=[0.1, 0.4, 0.7])

    simulation.run_equilibrium_time_simulation(measure="Mean", sample_size=100, number_of_steps=1000)

    data = simulation.measurements_to_dataframe()
    simulation.compute_equilibrium_time(data=data, sample_size=100, number_of_steps=1000, eval_confidence_range=0.1,
                                        eval_confidence_window=10, measure="Mean", fma=fma)

    # simulation.run_expectation_value_simulation(
    #     measures=["Mean", "Config"], n_measurements=1000, n_steps_equilibrium=100, n_steps_autocorrelation=100)
    #
    # data = simulation.measurements_to_dataframe()
    # simulation.compute_expectation_value(measures=["Mean", "Config"],
    #     data=data, eval_n_means_bootstrap=0)  # measures=["Mean", "SecondMoment"])
    print(1)
    # simulation.load_config_data()
