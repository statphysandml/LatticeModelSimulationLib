#ifndef LATTICEMODELSIMULATIONLIB_LINK_SAMPLER_HPP
#define LATTICEMODELSIMULATIONLIB_LINK_SAMPLER_HPP


#include "sampler_base.hpp"


namespace lm_impl {
    namespace link_lattice_system {
        struct LinkSampler : public lm_impl::sampler::SamplerBase<LinkSampler>
        {
            explicit LinkSampler(json params):
                SamplerBase<LinkSampler>(params),
                eps_(get_entry<double>("eps", 0.95))
            {}

            LinkSampler(const double eps=0.95):
                LinkSampler(json{{"eps", eps}})
            {}

            template<typename T>
            T cold_sample() {
                return T("identity");
            }

            template<typename T>
            T random_sample() {
                return T("random");
            }

            template<typename T>
            T propose_sample(T link) {
                return link * T(eps_);
            }

            double get_eps() const
            {
                return eps_;
            }

            double eps_;
        };
    }
}

#endif