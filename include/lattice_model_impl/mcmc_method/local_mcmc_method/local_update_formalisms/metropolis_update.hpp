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
                MCMCMethodBase<MetropolisUpdate<Model, Sampler>>(params),
                eps_(this->template get_entry<double>("eps", 0.0)),
                sampler_(Sampler(eps_))
            {
                rand_ = std::uniform_real_distribution<double>(0, 1);
            }

            explicit MetropolisUpdate(const double eps=0.0):
                MetropolisUpdate(json{{"eps", eps}})
            {}

            template<typename System>
            void initialize(const System &system) {
                model_ptr_ = &system.get_mcmc_model();
            }

            template<typename T>
            T random_state() {
                return sampler_.template random_state<T>();
            };

            template<typename T>
            T cold_state() {
                return sampler_.template cold_state<T>();
            };

            template<typename T>
            T operator()(const T site) {
                T proposed_site = sampler_.propose_state(site);
                if (rand_(mcmc::util::g_gen) <
                    std::min(1.0, std::exp(-1.0 * (model_ptr_->get_potential(proposed_site) - model_ptr_->get_potential(site)))))
                    return proposed_site;
                else
                    return site;
            }

            template<typename T>
            T operator()(const T site, const std::vector<T*> neighbours) {
                T proposed_site = sampler_.propose_state(site);
                if (rand_(mcmc::util::g_gen) < std::min(1.0, std::exp(-1.0 *
                                                                   (model_ptr_->get_potential(proposed_site, neighbours) -
                                                                    model_ptr_->get_potential(site, neighbours)))))
                    return proposed_site;
                else
                    return site;
            }

        protected:
            const Model* model_ptr_;

            Sampler sampler_;
            std::uniform_real_distribution<double> rand_;
            double eps_;
        };

    }
}

#endif //MAIN_METROPOLIS_UPDATE_HPP
