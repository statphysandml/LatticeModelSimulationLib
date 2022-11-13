#ifndef LATTICEMODELIMPLEMENTATIONS_DUMMY_UPDATE_DYNAMICS_HPP
#define LATTICEMODELIMPLEMENTATIONS_DUMMY_UPDATE_DYNAMICS_HPP


#include <lattice_model_impl/update_dynamics/update_dynamics_base.hpp>


namespace lm_impl {
    namespace update_dynamics {

        struct DummyUpdateDynamics : public UpdateDynamicsBase<DummyUpdateDynamics>
        {
            explicit DummyUpdateDynamics(const json params):
                UpdateDynamicsBase(params) {}

            explicit DummyUpdateDynamics():
                DummyUpdateDynamics(json{}) {}

            /* static std::string name() {
                return "DummyUpdateDynamics";
            } */

            template<typename Site>
            void initialize_update(const Site &site) {}

            template<typename Site>
            void update(Site &site, uint measure_interval = 1) {}
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_DUMMY_UPDATE_DYNAMICS_HPP
