#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_UPDATE_ON_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_UPDATE_ON_HPP


#include "langevin_update_base.hpp"


namespace lm_impl {
    namespace mcmc_method {
        template<typename Model>
        class ComplexLangevinONModelUpdate : public LangevinUpdateBase<ComplexLangevinONModelUpdate<Model>, Model> {
        public:
            using LangevinUpdateBase<ComplexLangevinONModelUpdate<Model>, Model>::LangevinUpdateBase;

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++) {
                    new_site(i).real(site(i).real() - epsilon * drift_term(i).real() + sqrt2epsilon * this->normal_(mcmc::util::g_gen));
                    new_site(i).imag(site(i).imag() - epsilon * drift_term(i).imag());
                }
                return this->model_ptr_->normalize_state(new_site);
            }
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_UPDATE_ON_HPP
