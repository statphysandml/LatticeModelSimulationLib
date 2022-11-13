#ifndef LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP


#include <lattice_model_impl/update_dynamics/update_dynamics_base.hpp>


namespace lm_impl {
    namespace update_dynamics {
        struct SiteSimpleUpdate : public UpdateDynamicsBase<SiteSimpleUpdate> {
            explicit SiteSimpleUpdate(const json params):
                UpdateDynamicsBase(params)
            {}

            explicit SiteSimpleUpdate() : SiteSimpleUpdate(json{}) {}

            static const std::string type() {
                return "SiteSimpleUpdate";
            }

            template<typename Site>
            void update(Site &site, uint measure_interval = 1) {
                for (uint k = 0; k < measure_interval; k++) {
                    site.get_system_representation() = update_system_site(*site.get_mcmc_method(), site.get_system_representation());
                }
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP
