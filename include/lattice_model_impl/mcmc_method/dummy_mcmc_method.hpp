#ifndef LATTICEMODELIMPLEMENTATIONS_DUMMY_mcmc_method_HPP
#define LATTICEMODELIMPLEMENTATIONS_DUMMY_mcmc_method_HPP


#include <lattice_model_impl/mcmc_method/mcmc_method_base.hpp>


namespace lm_impl {
    namespace mcmc_method {
        class DummyMCMCMethod : public MCMCMethodBase<DummyMCMCMethod> {
        public:
            template<typename Model>
            explicit DummyMCMCMethod(json params):
                MCMCMethodBase<DummyMCMCMethod>(params)
            {}

            template<typename T>
            T operator()(const T site) {
                return site;
            }

            template<typename T>
            T operator()(const T site, const std::vector<T*> neighbours) {
                return site;
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_DUMMY_mcmc_method_HPP
