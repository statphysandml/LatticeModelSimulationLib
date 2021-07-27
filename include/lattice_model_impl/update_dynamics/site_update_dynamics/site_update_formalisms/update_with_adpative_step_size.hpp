//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct UpdateWithAdpativeStepsize;

        struct UpdateWithAdpativeStepsizeParameters : UpdateDynamicsBaseParameters {
            explicit UpdateWithAdpativeStepsizeParameters(const json params_) : UpdateDynamicsBaseParameters(params_) {
                thermalization_langevin_time_interval = get_entry<double>("thermalization_langevin_time_interval", 10.0);
                langevin_time_measure_interval = get_entry<double>("langevin_time_measure_interval", 1.0);
            }

            explicit UpdateWithAdpativeStepsizeParameters(double thermalization_langevin_time_interval_, double langevin_time_measure_interval_) :
                    UpdateWithAdpativeStepsizeParameters(json {
                        {"thermalization_langevin_time_interval", thermalization_langevin_time_interval_},
                        {"langevin_time_measure_interval", langevin_time_measure_interval_}})
            {}

            static std::string name() {
                return "UpdateWithAdpativeStepsize";
            }

            typedef UpdateWithAdpativeStepsize UpdateDynamics;

            double thermalization_langevin_time_interval;
            double langevin_time_measure_interval;
        };

        struct UpdateWithAdpativeStepsize : public UpdateDynamicsBase<UpdateWithAdpativeStepsize> {
            explicit UpdateWithAdpativeStepsize(const UpdateWithAdpativeStepsizeParameters &sp_) : sp(sp_) {}

            template<typename Site>
            void initialize_update(const Site &site) {
                thermalized = false;
            }

            template<typename Site>
            void update(Site &lattice, uint measure_interval = 1) {
                if (thermalized)
                    parallel_update_with_adpative_stepsize(lattice, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(lattice);
            }

            template<typename Site>
            void thermalization_phase_with_adpative_stepsize(Site &site) {
                double max_epsilon = site.get_update_formalism().get_stepsize();
                int n_thermalization_steps = int(sp.thermalization_langevin_time_interval / max_epsilon);
                std::cout << "Perform " << n_thermalization_steps << " thermalization steps with a step width of " << max_epsilon << std::endl;
                for (int k = 0; k < n_thermalization_steps; k++) {
                    KExpectation += std::fabs(site.get_update_formalism().estimate_drift_term(site.get_system_representation()));
                    site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                           site.get_system_representation());
                }
                KExpectation /= n_thermalization_steps;
                thermalized = true;
            }


            template<typename Site>
            void parallel_update_with_adpative_stepsize(Site &site, uint measure_interval = 1) {
                double max_epsilon = site.get_update_formalism().get_stepsize();

                // Loop for evolving measure_interval * sp.langevin_time_measure_interval in the Langevin time
                for (uint k = 0; k < measure_interval; k++) {
                    uint break_out_counter = 0;
                    while(langevin_time < sp.langevin_time_measure_interval and break_out_counter < 10000000) {
                        const double KMax = std::fabs(site.get_update_formalism().estimate_drift_term(site.get_system_representation()));
                        double epsilon = std::min(max_epsilon, max_epsilon * KExpectation / KMax);
                        if(epsilon > 1.0) {
                            std::cout << "Detected epsilon > 1 in adaptive stepsize" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        site.get_system_representation() = update_lattice_site(site.get_update_formalism(),
                                                                               site.get_system_representation(), epsilon);
                        langevin_time += epsilon;
                        break_out_counter += 1;
                        // std::cout << "Epsilon " << epsilon << std::endl;
                    }
                    langevin_time -= sp.langevin_time_measure_interval;

                    // std::cout << "Langevin time" << langevin_time << std::endl;

                    if(break_out_counter == 10000000)
                    {
                        std::cout << "Too small stepsize, break_out_counter == 1000000 for one Langevin time step" << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                }


            }

            const UpdateWithAdpativeStepsizeParameters &sp;

            // For adaptive step size
            double KExpectation = 0;
            bool thermalized = 0;

            double langevin_time = 0;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
