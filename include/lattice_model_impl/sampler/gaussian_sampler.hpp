#ifndef MAIN_GAUSSIAN_SAMPLER_HPP
#define MAIN_GAUSSIAN_SAMPLER_HPP


#include <mcmc_simulation/util/random.hpp>

#include "sampler_base.hpp"


namespace lm_impl {
    namespace sampler {

        /** @brief Gaussian distribution sampler with variance 2 * eps
         *
         * Sampler function for sampling, evaluating and integrating a Gaussian distribution.
         */
        struct GaussianSampler : public SamplerBase<GaussianSampler> {
            
            explicit GaussianSampler(json params):
                SamplerBase<GaussianSampler>(params),
                eps_(this->template get_entry<double>("eps", 0.1))
            {
                normal_ = std::normal_distribution<double>(0, 1);
            }

            explicit GaussianSampler(const double eps=0.1):
                GaussianSampler(json{{"eps", eps}})
            {}

            template<typename T>
            T random_sample() {
                return std::sqrt(2 * eps_) * normal_(mcmc::util::g_gen);
            }

            template<typename T>
            T cold_sample() {
                return T(0);
            }

            template<typename T>
            T proposal_sample(T site) {
                return site + std::sqrt(2 * eps_) * normal_(mcmc::util::g_gen);
            }

            double get_eps() const {
                return eps_;
            }

            template<typename T>
            T evaluate(T site) const {
                return 1.0 / sqrt(2 * M_PI) * std::exp(-1.0 * std::pow(site, 2) / 2.0);
            }

            template<typename T>
            std::pair<double, double> get_integration_bounds(const T &site) const {
                return std::pair<double, double>(-1.0, 1.0);
            }

            struct transformer_func {
                double operator()(const double val) {
                    return std::log((1.0 + val) / (1.0 - val));
                }
            };

            double jacobian(const double x) {
                return -2.0 / (std::pow(x, 2.0) - 1.0);
            }

            double eps_;
            std::normal_distribution<double> normal_;

            transformer_func transformer_;
        };
    }
}

#endif //MAIN_GAUSSIAN_SAMPLER_HPP
