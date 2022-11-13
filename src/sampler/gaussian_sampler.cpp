#include <lattice_model_impl/sampler/gaussian_sampler.hpp>

namespace lm_impl {
    namespace sampler {

        template<>
        std::complex<double> GaussianSampler::random_sample() {
            return {std::sqrt(2 * eps_) * normal_(mcmc::util::random::g_gen), 0.1};
        }

        template<>
        std::complex<double> GaussianSampler::propose_sample(std::complex<double> site) {
            return site + std::complex<double>{std::sqrt(2 * eps_) * normal_(mcmc::util::random::g_gen), 0.0};
        }

    }
}