//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct ParallelUpdateWithAdpativeStepsize;

        struct ParallelUpdateWithAdpativeStepsizeParameters : UpdateDynamicsBaseParameters {
            explicit ParallelUpdateWithAdpativeStepsizeParameters(const json params_) : UpdateDynamicsBaseParameters(
                    params_) {
                thermalization_langevin_time_interval = get_entry<double>("thermalization_langevin_time_interval", 10.0);
                langevin_time_measure_interval = get_entry<double>("langevin_time_measure_interval", 1.0);
            }

            explicit ParallelUpdateWithAdpativeStepsizeParameters(double thermalization_langevin_time_interval_, double langevin_time_measure_interval_) :
                    ParallelUpdateWithAdpativeStepsizeParameters(json {
                            {"thermalization_langevin_time_interval", thermalization_langevin_time_interval_},
                            {"langevin_time_measure_interval", langevin_time_measure_interval_}})
            {}

            static std::string name() {
                return "ParallelUpdateWithAdpativeStepsize";
            }

            typedef ParallelUpdateWithAdpativeStepsize UpdateDynamics;

            double thermalization_langevin_time_interval;
            double langevin_time_measure_interval;
        };

        struct ParallelUpdateWithAdpativeStepsize : public UpdateDynamicsBase<ParallelUpdateWithAdpativeStepsize> {
            explicit ParallelUpdateWithAdpativeStepsize(const ParallelUpdateWithAdpativeStepsizeParameters &lp_) : lp(
                    lp_) {}

            template<typename Lattice>
            void initialize_update(const Lattice &lattice) {
                thermalized = false;
            }

            template<typename Lattice>
            void update(Lattice &lattice, uint measure_interval = 1) {
                if (thermalized)
                    parallel_update_with_adpative_stepsize(lattice, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(lattice);
            }

            template<typename Lattice>
            void thermalization_phase_with_adpative_stepsize(Lattice &lattice, uint measure_interval = 1) {
                // static_assert(detail::is_updateable<T, typename UpdateFormalismParameters::MCMCMethod>::value, "is not estimate_drift_term");

                double max_epsilon = lattice.get_update_formalism().get_stepsize();
                int n_thermalization_steps = int(lp.thermalization_langevin_time_interval / max_epsilon);
                std::cout << "Perform " << n_thermalization_steps << " thermalization steps with a step width of " << max_epsilon << std::endl;

                for (int k = 0; k < n_thermalization_steps; k++) {

                    double KMax = 0;
                    for (uint i = 0; i < lattice.size(); i++) {
                        const double K = std::fabs(
                                lattice.get_update_formalism().estimate_drift_term(lattice[i],
                                                                                   lattice.neighbours_at(i)));
                        if (K > KMax)
                            KMax = K;
                    }
                    KExpectation += KMax;

                    std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.size(),
                                                                             typename Lattice::SiteType(0));

                    // #pragma omp parallel for
                    for (uint i = 0; i < lattice.size(); i++) {
                        lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i],
                                                                  lattice.neighbours_at(i));
                    }

                    // Rewrite?
                    auto &lattice_grid = lattice.get_system_representation();
                    lattice_grid = lattice_grid_new;
                }

                KExpectation /= n_thermalization_steps;
                thermalized = true;
            }


            template<typename Lattice>
            void parallel_update_with_adpative_stepsize(Lattice &lattice, uint measure_interval = 1) {
                double max_epsilon = lattice.get_update_formalism().get_stepsize();

                for (uint j = 0; j < measure_interval; j++) {
                    uint break_out_counter = 0;
                    while(langevin_time < lp.langevin_time_measure_interval and break_out_counter < 10000000) {
                        double KMax = 0;
                        for (uint i = 0; i < lattice.size(); i++) {
                            const double K = std::fabs(
                                    lattice.get_update_formalism().estimate_drift_term(lattice[i],
                                                                                       lattice.neighbours_at(i)));
                            if (K > KMax)
                                KMax = K;
                        }

                        double epsilon = std::min(max_epsilon, max_epsilon * KExpectation / KMax);

                        std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.size(),
                                                                                 typename Lattice::SiteType(0));

                        // #pragma omp parallel for
                        for (uint i = 0; i < lattice.size(); i++) {
                            lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i],
                                                                      lattice.neighbours_at(i), epsilon);
                        }

                        // Rewrite?
                        auto &lattice_grid = lattice.get_system_representation();
                        lattice_grid = lattice_grid_new;

                        langevin_time += epsilon;
                        break_out_counter += 1;
                    }

                    langevin_time -= lp.langevin_time_measure_interval;

                    if(break_out_counter == 10000000)
                    {
                        std::cerr << "Too small stepsize, break_out_counter == 1000000 for one Langevin time step" << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                }
            }

            const ParallelUpdateWithAdpativeStepsizeParameters &lp;

            // For adaptive step size
            double KExpectation = 0;
            bool thermalized = 0;

            double langevin_time = 0;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
