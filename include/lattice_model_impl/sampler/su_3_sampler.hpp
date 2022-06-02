#ifndef LATTICEMODELSIMULATIONLIB_SU3_SAMPLER_HPP
#define LATTICEMODELSIMULATIONLIB_SU3_SAMPLER_HPP


#include "sampler_base.hpp"
#include <math.h>




namespace lm_impl {
    namespace sampler {
        struct SU3Sampler : public lm_impl::sampler::SamplerBase<SU3Sampler>
        {
            explicit SU3Sampler(json params):
                SamplerBase<SU3Sampler>(params),
                eps_(get_entry<double>("eps", 0.1))
            {
                normal_ = std::normal_distribution<double>(0, 1);
            }

            SU3Sampler(const double eps=0.1):
                SU3Sampler(json{{"eps", eps}})
            {}            

            template<typename T>
            T random_sample() {
                double pi = M_PI;
                std::uniform_real_distribution<double> dis(-pi, pi);


                double phi1 = dis(mcmc::util::g_gen);
                double phi2 = dis(mcmc::util::g_gen);

                T new_link(std::exp(1i*phi1), 0, 0,
                            0, std::exp(1i*phi2), 0,
                            0, 0, std::exp(-1i*(phi1 + phi2)));
                
                return new_link;
            }

            template<typename T>
            T cold_sample() {
                return T(1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0);
            }

            //template<typename T>
            //T propose_sample(T site) {
            //    return site + random_sample<T>();
            //}

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