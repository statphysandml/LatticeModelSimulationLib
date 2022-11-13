#ifndef MAIN_LANGEVIN_UPDATE_HPP
#define MAIN_LANGEVIN_UPDATE_HPP


#include <lattice_model_impl/mcmc_method/local_mcmc_method/local_update_formalisms/langevin_update_base.hpp>


namespace lm_impl {
    namespace mcmc_method {
        template<typename Model>
        class LangevinUpdate : public LangevinUpdateBase<LangevinUpdate<Model>, Model> {
        public:
            using LangevinUpdateBase<LangevinUpdate<Model>, Model>::LangevinUpdateBase;

            static const std::string type() {
                return "LangevinUpdate";
            }

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                return this->model_ptr_->normalize_state(site - epsilon * drift_term + sqrt2epsilon * this->normal_(mcmc::util::random::g_gen));
            }
        };
    }
}

#endif //MAIN_LANGEVIN_UPDATE_HPP
