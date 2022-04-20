#include "../../../include/lattice_model_impl/sampler/hat_function_sampler.hpp"

namespace lm_impl {
    namespace sampler {

        template<>
        std::complex<double> HatFunctionSampler::random_state() {
            auto r = uniform(mcmc::util::g_gen);
            if (r < 0.5)
                return {-1.0 * eps * (1.0 - sqrt(2 * r)), -0.6};
            else
                return {1.0 * eps * (1.0 - sqrt(2.0 * (1.0 - r))), -0.6};
        }

        template<>
        std::complex<double> HatFunctionSampler::propose_state(std::complex<double> site) {
            auto r = uniform(mcmc::util::g_gen);
            if (r < 0.5)
                return site + std::complex<double>{-1.0 * eps * (1.0 - sqrt(2 * r)), 0};
            else
                return site + std::complex<double>{1.0 * eps * (1.0 - sqrt(2.0 * (1.0 - r))), 0};
        }

    }
}