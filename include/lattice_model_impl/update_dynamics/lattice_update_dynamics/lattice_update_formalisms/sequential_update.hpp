#ifndef LATTICEMODELIMPLEMENTATIONS_SEQUENTIAL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SEQUENTIAL_UPDATE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {
        struct SequentialUpdate : public UpdateDynamicsBase<SequentialUpdate> {
            explicit SequentialUpdate(const json params):
                UpdateDynamicsBase(params)
            {}

            explicit SequentialUpdate() : SequentialUpdate(json{}) {}

            template<typename System>
            void initialize(System &system) {
                uniint_ = std::uniform_int_distribution<int>(0, system.size() - 1);
            }

            template<typename System>
            void update(System &system, uint measure_interval = 1) {
                for (size_t k = 0; k < measure_interval; k++) {
                    for (uint j = 0; j < system.size(); j++) {
                        int i = uniint_(mcmc::util::g_gen);
                        system[i] = update_system_site(system.get_mcmc_method(), system[i], system.neighbours_at(i));
                        // const double K = std::fabs(update_formalism->estimate_drift_term(system[i], system.neighbours_at[i]));
                    }
                }
            }

            std::uniform_int_distribution<int> uniint_;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SEQUENTIAL_UPDATE_HPP
