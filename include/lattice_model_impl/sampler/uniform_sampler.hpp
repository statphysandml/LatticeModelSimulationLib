#ifndef MAIN_UNIFORM_SAMPLER_HPP
#define MAIN_UNIFORM_SAMPLER_HPP


#include <mcmc/mcmc_simulation/util/random.hpp>

#include <lattice_model_impl/sampler/sampler_base.hpp>


namespace lm_impl {
    namespace sampler {

        /** @brief Uniform distribution sampler with witdh 2 * eps
         *
         * Sampler function for sampling, evaluating and integrating a uniform distribution.
         */
        struct UniformSampler : public lm_impl::sampler::SamplerBase<UniformSampler>
        {
            explicit UniformSampler(json params):
                lm_impl::sampler::SamplerBase<UniformSampler>(params),
                eps_(this->template get_entry<double>("eps", 0.1))
            {
                uniform_ = std::uniform_real_distribution<double>(-1.0, 1.0);
            }


            UniformSampler(const double eps=0.1):
                UniformSampler(json{{"eps", eps}})
            {}
            
            static const std::string type() {
                return "UniformSampler";
            }

            template<typename T>
            T random_sample() {
                return eps_ * uniform_(mcmc::util::random::g_gen);
            }

            template<typename T>
            T cold_sample() {
                return T(0);
            }

            template<typename T>
            T propose_sample(T site) {
                return site + eps_ * uniform_(mcmc::util::random::g_gen);
            }

            double get_eps() const {
                return eps_;
            }

            template<typename T>
            std::pair<double, double> get_integration_bounds(const T &site) const {
                return std::pair<double, double>(site.real() - eps_, site.real() + eps_);
            }

            struct transformer_func {
                double operator()(const double val) {
                    return val;
                }
            };

            double jacobian(const double x) {
                return 1.0;
            }

            double eps_;
            std::uniform_real_distribution<double> uniform_;

            transformer_func transformer_;
        };
    }
}

#endif //MAIN_UNIFORM_SAMPLER_HPP
