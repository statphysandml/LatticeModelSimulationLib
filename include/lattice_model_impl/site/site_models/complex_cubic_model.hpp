#ifndef MAIN_COMPLEX_CUBIC_MODEL_HPP
#define MAIN_COMPLEX_CUBIC_MODEL_HPP

#include <lattice_model_impl/lattice/mcmc_model_base.hpp>
#include <mcmc/mcmc_simulation/util/random.hpp>


namespace lm_impl {
    namespace site_system {
        class ComplexCubicModel : public lm_impl::model::MCMCModelBase<ComplexCubicModel> {
        public:
            explicit ComplexCubicModel(const json params):
                MCMCModelBase(params)
            {}

            explicit ComplexCubicModel() : ComplexCubicModel(json{}) {}

            static const std::string type() {
                return "ComplexCubicModel";
            }

            static std::complex<double> get_drift_term(const std::complex<double> site) {
                return {-2.0 * site.real() * site.imag(), -1.0 * (std::pow(site.imag(), 2) - std::pow(site.real(), 2))};
            }

            static std::complex<double> get_second_order_drift_term(const std::complex<double> site) {
                return {-2.0 * site.imag(), 2.0 * site.real()};
            }

            static std::complex<double> get_potential(const std::complex<double> site) {
                // return std::complex<double>{0, 1} * std::pow(site, 3) / 3.0;
                return {-1.0 * std::pow(site.real(), 2) * site.imag() + std::pow(site.imag(), 3) / 3.0,
                        std::pow(site.real(), 3) / 3.0 - std::pow(site.imag(), 2) * site.real()};
            }
        };
    }
}

#endif //MAIN_COMPLEX_CUBIC_MODEL_HPP
