#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {
        struct ParallelUpdateWithAdpativeStepsize : public UpdateDynamicsBase<ParallelUpdateWithAdpativeStepsize> {
            explicit ParallelUpdateWithAdpativeStepsize(const json params):
                UpdateDynamicsBase(params),
                thermalization_langevin_time_interval_(get_entry<double>("thermalization_langevin_time_interval", 10.0)),
                langevin_time_measure_interval_(get_entry<double>("langevin_time_measure_interval", 1.0))
            {}

            explicit ParallelUpdateWithAdpativeStepsize(double thermalization_langevin_time_interval=10.0, double langevin_time_measure_interval=1.0):
                ParallelUpdateWithAdpativeStepsize(json {
                    {"thermalization_langevin_time_interval", thermalization_langevin_time_interval},
                    {"langevin_time_measure_interval", langevin_time_measure_interval}})
            {}

            template<typename System>
            void initialize_update(System &system) {
                thermalized_ = false;
            }

            template<typename System>
            void update(System &system, uint measure_interval = 1) {
                if (thermalized_)
                    parallel_update_with_adpative_stepsize(system, measure_interval);
                else
                    thermalization_phase_with_adpative_stepsize(system);
            }

            template<typename System>
            void thermalization_phase_with_adpative_stepsize(System &system) {
                double max_epsilon = system.get_mcmc_method().get_stepsize();
                int n_thermalization_steps = int(thermalization_langevin_time_interval_ / max_epsilon);
                std::cout << "Perform " << n_thermalization_steps << " thermalization steps with a step width of " << max_epsilon << std::endl;

                for (int k = 0; k < n_thermalization_steps; k++) {

                    double KMax = 0;
                    for (uint i = 0; i < system.size(); i++) {
                        const double K = std::fabs(
                                system.get_mcmc_method().estimate_drift_term(system[i], system.neighbours_at(i)));
                        if (K > KMax)
                            KMax = K;
                    }
                    KExpectation_ += KMax;

                    std::vector<typename System::SiteType> system_grid_new(system.size(), typename System::SiteType(0));

                    for (uint i = 0; i < system.size(); i++) {
                        system_grid_new[i] = update_system_site(system.get_mcmc_method(), system[i],
                                                                  system.neighbours_at(i));
                    }

                    auto &system_grid = system.get_system_representation();
                    system_grid = system_grid_new;
                }

                KExpectation_ /= n_thermalization_steps;
                thermalized_ = true;
            }

            template<typename System>
            void parallel_update_with_adpative_stepsize(System &system, uint measure_interval = 1) {
                double max_epsilon = system.get_mcmc_method().get_stepsize();

                for (uint j = 0; j < measure_interval; j++) {
                    uint break_out_counter = 0;
                    while(langevin_time_ < langevin_time_measure_interval_ and break_out_counter < 10000000) {
                        double KMax = 0;
                        for (uint i = 0; i < system.size(); i++) {
                            const double K = std::fabs(
                                    system.get_mcmc_method().estimate_drift_term(system[i], system.neighbours_at(i)));
                            if (K > KMax)
                                KMax = K;
                        }

                        double epsilon = std::min(max_epsilon, max_epsilon * KExpectation_ / KMax);

                        std::vector<typename System::SiteType> system_grid_new(system.size(), typename System::SiteType(0));

                        for (uint i = 0; i < system.size(); i++) {
                            system_grid_new[i] = update_system_site(system.get_mcmc_method(), system[i],
                                                                      system.neighbours_at(i), epsilon);
                        }

                        auto &system_grid = system.get_system_representation();
                        system_grid = system_grid_new;

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

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
