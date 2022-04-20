//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_METROPOLIS_UPDATE_HPP
#define MAIN_METROPOLIS_UPDATE_HPP


#include "../../mcmc_method_base.hpp"


namespace lm_impl {
    namespace mcmc_method {

        template<typename Model, typename Sampler>
        class MetropolisUpdate : public MCMCMethodBase<MetropolisUpdate<Model, Sampler>> {
        public:
            explicit MetropolisUpdate(json params) :
                MCMCMethodBase<MetropolisUpdate<Model, Sampler>>(params)
            {
                rand_ = std::uniform_real_distribution<double>(0, 1);
            }

            explicit MetropolisUpdate():
                MetropolisUpdate(json{})
            {}

            template<typename System>
            void initialize(System &system) {
                sampler_ptr_ = system.get_sampler().get();
                model_ptr_ = system.get_mcmc_model().get();
            }

            template<typename T>
            T operator()(const T site) {
                T proposed_site = sampler_ptr_->template propose_sample(site);
                if(rand_(mcmc::util::g_gen) <
                    std::min(1.0, std::exp(-1.0 * (model_ptr_->get_potential(proposed_site) - model_ptr_->get_potential(site)))))
                    return this->model_ptr_->normalize_state(proposed_site);
                else
                    return site;
            }

            template<typename T>
            T operator()(const T site, const std::vector<T*> neighbours) {
                T proposed_site = sampler_ptr_->template propose_sample(site);
                if(rand_(mcmc::util::g_gen) < std::min(
                        1.0, std::exp(
                            -1.0 * (model_ptr_->get_potential(proposed_site, neighbours) -
                                    model_ptr_->get_potential(site, neighbours))
                        )
                    )
                )
                    return this->model_ptr_->normalize_state(proposed_site);
                else
                    return site;
            }

        private:
            Sampler* sampler_ptr_;
            Model* model_ptr_;
            
            std::uniform_real_distribution<double> rand_;
        };

    }
}

#endif //MAIN_METROPOLIS_UPDATE_HPP
