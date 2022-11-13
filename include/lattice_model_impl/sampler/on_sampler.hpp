#ifndef LATTICEMODELSIMULATIONLIB_ON_SAMPLER_HPP
#define LATTICEMODELSIMULATIONLIB_ON_SAMPLER_HPP


#include <lattice_model_impl/sampler/sampler_base.hpp>


namespace lm_impl {
    namespace lattice_system {
        struct ONSampler : public lm_impl::sampler::SamplerBase<ONSampler>
        {
            explicit ONSampler(json params):
                SamplerBase<ONSampler>(params),
                eps_(get_entry<double>("eps", 0.1))
            {
                normal_ = std::normal_distribution<double>(0, 1);
            }

            explicit ONSampler(const double eps=0.1):
                ONSampler(json{{"eps", eps}})
            {}            

            static const std::string type() {
                return "ONSampler";
            }

            template<typename T>
            T random_sample() {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++)
                    new_site(i) += std::sqrt(2 * eps_) * normal_(mcmc::util::random::g_gen);
                return new_site;
            }

            template<typename T>
            T cold_sample() {
                return T(0);
            }

            template<typename T>
            T propose_sample(T site) {
                return site + random_sample<T>();
            }

            double get_eps() const
            {
                return eps_;
            }

            double eps_;
            std::normal_distribution<double> normal_;
        };
    }
}

#endif