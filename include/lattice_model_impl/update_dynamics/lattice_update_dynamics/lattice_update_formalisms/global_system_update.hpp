//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP


#include "../../update_dynamics_base.hpp"


namespace lm_impl {
    namespace update_dynamics {

        struct GlobalSystemUpdate : public UpdateDynamicsBase<GlobalSystemUpdate> {
            explicit GlobalSystemUpdate(const json params) : UpdateDynamicsBase(params) {}

            explicit GlobalSystemUpdate() : GlobalSystemUpdate(json{}) {}

            template<typename System>
            void update(System &system, uint measure_interval = 1) {
                for (uint k = 0; k < measure_interval; k++) {
                    global_system_update(system.get_mcmc_method(), system);
                }
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP
