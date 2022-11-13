#include <lattice_model_impl/sampler/uniform_sampler.hpp>

namespace lm_impl {
    namespace sampler {

        template<>
        std::complex<double> UniformSampler::random_sample() {
            return {eps_ * uniform_(mcmc::util::random::g_gen), 0.1};
        }

        template<>
        std::complex<double> UniformSampler::propose_sample(std::complex<double> site) {
            return site + std::complex<double>{eps_ * uniform_(mcmc::util::random::g_gen), 0.0};
        }

    }
}