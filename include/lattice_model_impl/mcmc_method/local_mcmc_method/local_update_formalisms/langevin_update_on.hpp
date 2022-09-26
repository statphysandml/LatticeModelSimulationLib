//
// Created by lukas on 20.03.21.
//

#ifndef LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_ON_HPP
#define LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_ON_HPP


#include "../../mcmc_method_base.hpp"
#include "../../../lattice/lattice_models/on_model.hpp"


namespace lm_impl {
    namespace mcmc_method {
        template<typename Model>
        class LangevinONModelUpdate : public LangevinUpdateBase<LangevinONModelUpdate<Model>, Model> {
        public:
            using LangevinUpdateBase<LangevinONModelUpdate<Model>, Model>::LangevinUpdateBase;

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++) {
                    new_site(i) = site(i) - epsilon * drift_term(i) + sqrt2epsilon * this->normal_(mcmc::util::g_gen);
                }
                return this->model_ptr_->normalize_state(new_site);
            }
        };
    }
}

#endif //LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_ON_HPP
