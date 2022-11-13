#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP


#include <lattice_model_impl/update_dynamics/update_dynamics_base.hpp>


namespace lm_impl {
    namespace update_dynamics {
        struct ParallelUpdate : public UpdateDynamicsBase<ParallelUpdate> {
            explicit ParallelUpdate(const json params):
                UpdateDynamicsBase(params)
            {}

            static const std::string type() {
                return "ParallelUpdate";
            }

            template<typename System>
            void update(System &system, uint measure_interval = 1) {
                for (uint j = 0; j < measure_interval; j++) {
                    std::vector<typename System::SiteType> system_grid_new(system.size(), typename System::SiteType(0));

                    for (uint i = 0; i < system.size(); i++) {
                        system_grid_new[i] = update_system_site(system.get_update_formalism(), system[i], system.neighbours_at(i));
                    }

                    auto &system_grid = system.get_system_representation();
                    system_grid = system_grid_new;
                }
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP
