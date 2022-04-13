#ifndef MAIN_COMPLEX_LANGEVIN_UPDATE_HPP
#define MAIN_COMPLEX_LANGEVIN_UPDATE_HPP


#include "langevin_update_base.hpp"


namespace lm_impl {
    namespace mcmc_method {
        template<typename Model>
        class ComplexLangevinUpdate : public LangevinUpdateBase<ComplexLangevinUpdate<Model>, Model> {
        public:
            using LangevinUpdateBase<ComplexLangevinUpdate<Model>, Model>::LangevinUpdateBase;

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                T new_site = {site.real(), site.imag()};
                new_site.real(site.real() - epsilon * drift_term.real() + sqrt2epsilon * this->normal_(mcmc::util::g_gen));
                new_site.imag(site.imag() - epsilon * drift_term.imag());

                std::cout << new_site.real() << " + i " <<  new_site.imag() << std::endl;

                return this->model_ptr_->normalize_state(new_site);
            }
        };
    }
}

#endif //MAIN_COMPLEX_LANGEVIN_UPDATE_HPP
