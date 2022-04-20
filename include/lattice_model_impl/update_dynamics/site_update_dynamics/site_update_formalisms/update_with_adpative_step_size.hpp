#ifndef LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {
        struct UpdateWithAdpativeStepsize : public UpdateDynamicsBase<UpdateWithAdpativeStepsize> {
            explicit UpdateWithAdpativeStepsize(const json params):
                UpdateDynamicsBase(params),
                thermalization_langevin_time_interval_(get_entry<double>("thermalization_langevin_time_interval", 100.0)),
                langevin_time_measure_interval_(get_entry<double>("langevin_time_measure_interval", 1.0))
            {}

            explicit UpdateWithAdpativeStepsize(double thermalization_langevin_time_interval=10.0, double langevin_time_measure_interval=1.0):
                UpdateWithAdpativeStepsize(json {
                    {"thermalization_langevin_time_interval", thermalization_langevin_time_interval},
                    {"langevin_time_measure_interval", langevin_time_measure_interval}})
            {}

            template<typename Site>
            void initialize(Site &site) {
                thermalized_ = false;
            }

            template<typename Site>
            void update(Site &site, uint measure_interval = 1) {
                if (thermalized_)
                    parallel_update_with_adpative_stepsize(site, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(site);
            }

            template<typename Site>
            void thermalization_phase_with_adpative_stepsize(Site &site) {
                double max_epsilon = site.get_mcmc_method()->get_stepsize();
                int n_thermalization_steps = int(thermalization_langevin_time_interval_ / max_epsilon);
                std::cout << "Perform " << n_thermalization_steps << " thermalization steps with a step width of " << max_epsilon << std::endl;
                for (int k = 0; k < n_thermalization_steps; k++) {
                    KExpectation_ += std::fabs(site.get_mcmc_method()->estimate_drift_term(site.get_system_representation()));
                    site.get_system_representation() = update_system_site(*site.get_mcmc_method(), site.get_system_representation());
                }
                KExpectation_ /= n_thermalization_steps;
                thermalized_ = true;
            }


            template<typename Site>
            void parallel_update_with_adpative_stepsize(Site &site, uint measure_interval = 1) {
                double max_epsilon = site.get_mcmc_method()->get_stepsize();

                // Loop for evolving measure_interval * langevin_time_measure_interval in the Langevin time
                for (uint k = 0; k < measure_interval; k++) {
                    uint break_out_counter = 0;
                    while(langevin_time_ < langevin_time_measure_interval_ and break_out_counter < 10000000) {
                        const double KMax = std::fabs(site.get_mcmc_method()->estimate_drift_term(site.get_system_representation()));
                        double epsilon = std::min(max_epsilon, max_epsilon * KExpectation_ / KMax);
                        if(epsilon > 1.0) {
                            std::cerr << "Detected epsilon > 1 in adaptive stepsize" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        site.get_system_representation() = update_system_site(*site.get_mcmc_method(), site.get_system_representation(), epsilon);
                        langevin_time_ += epsilon;
                        break_out_counter += 1;
                    }
                    langevin_time_ -= langevin_time_measure_interval_;

                    if(break_out_counter == 10000000)
                    {
                        std::cerr << "Too small stepsize, break_out_counter == 1000000 for one Langevin time step" << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                }
            }

            // For adaptive step size
            double KExpectation_ = 0;
            bool thermalized_ = 0;

            double langevin_time_ = 0;
            double thermalization_langevin_time_interval_;
            double langevin_time_measure_interval_;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_UPDATE_WITH_ADPATIVE_STEP_SIZE_HPP
