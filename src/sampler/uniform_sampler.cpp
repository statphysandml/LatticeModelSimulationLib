#include "../../../include/lattice_model_impl/sampler/uniform_sampler.hpp"

namespace lm_impl {
    namespace sampler {

        template<>
        std::complex<double> UniformSampler::random_state() {
            return {eps * uniform(mcmc::util::g_gen), 0.1};
        }

        template<>
        std::complex<double> UniformSampler::propose_state(std::complex<double> site) {
            return site + std::complex<double>{eps * uniform(mcmc::util::g_gen), 0.0};
        }

    }
}