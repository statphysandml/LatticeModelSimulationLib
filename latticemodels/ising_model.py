
class IsingModel:
    def __init__(self, beta, J, h, dimensions, lattice_update_formalism="Sequential", **kwargs):
        print("Calling Ising Model Base Head")
        super().__init__(**kwargs)
        print("Calling Ising Model Base Body")

        self.beta = beta
        self.J = J
        self.h = h
        self.dimensions = dimensions
        self.lattice_update_formalism = lattice_update_formalism

        self.measures = []
        self.parameters = None

    """ Possibility to include if conditions! """

    def initialize_parameters(self):
        from latticemodelimpl import IsingModelParameters
        self.model = IsingModelParameters(beta=self.beta, J=self.J, h=self.h)

        from latticemodelimpl import IsingModelMetropolisParameters
        self.mcmc_method = IsingModelMetropolisParameters()

        # Possibility to have this as a factory method!
        if self.lattice_update_formalism == "Sequential":
            from latticemodelimpl import SequentialUpdateParameters
            lattice_update_formalism_ = SequentialUpdateParameters()
        else:
            assert False, "Unkown lattice update formalism"

        from latticemodelimpl import IsingModelMetropolisLatticeParameters as LatticeParameters
        self.parameters = LatticeParameters.generate_parameters(
            self.model,
            self.mcmc_method,
            lattice_update_formalism_,
            self.dimensions,
            "nearest_neighbour",
            self.measures
        )

    def set_measures(self, measures):
        self.measures = measures

    def generate_mcmc_system(self):
        from latticemodelimpl import IsingModelMetropolisLattice as Lattice
        return Lattice(self.parameters)

    def _equilibrium_time_simulation_classes(self):
        from latticemodelimpl import IsingModelMetropolisEquilibriumTimeParameters as EquiTimeSimParams
        from latticemodelimpl import IsingModelMetropolisEquilibriumTime as EquiTimeSim
        return EquiTimeSimParams, EquiTimeSim

    def _correlation_time_simulation_classes(self):
        from latticemodelimpl import IsingModelMetropolisCorrelationTimeParameters as CorrTimeSimParams
        from latticemodelimpl import IsingModelMetropolisCorrelationTime as CorrTimeSim
        return CorrTimeSimParams, CorrTimeSim

    def _expectation_value_simulation_classes(self):
        from latticemodelimpl import IsingModelMetropolisExpectationValueParameters as ExpValSimParams
        from latticemodelimpl import IsingModelMetropolisExpectationValue as ExpValSim
        return ExpValSimParams, ExpValSim
