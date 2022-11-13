#ifndef LATTICEMODELIMPLEMENTATIONS_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_POLYNOMIAL_MODEL_HPP

#include <lattice_model_impl/lattice/mcmc_model_base.hpp>
#include <mcmc/mcmc_simulation/util/random.hpp>


namespace lm_impl {
    namespace site_system {
        class PolynomialModel : public lm_impl::model::MCMCModelBase<PolynomialModel> {
        public:
            explicit PolynomialModel(const json params):
                MCMCModelBase(params),
                sigma_(get_entry<double>("sigma", 1.0)),
                lambda_(get_entry<double>("lambda", 1.0)),
                h_(get_entry<double>("h", 0.0)) {}

            explicit PolynomialModel(double lambda=1.0, double sigma=1.0, double h=0.0):
                PolynomialModel(
                    json{
                            {"lambda", lambda},
                            {"sigma",  sigma},
                            {"h",      h}
                    }) {}

            static const std::string type() {
                return "PolynomialModel";
            }

            double get_potential(const double site) const {
                return 0.5 * sigma_ * std::pow(site, 2) + 0.25 * lambda_ * std::pow(site, 4) + h_ * site;
            }

            double get_drift_term(const double site) const {
                return sigma_ * site + lambda_ * std::pow(site, 3) + h_;
            }

        private:
            double sigma_;
            double lambda_;
            double h_;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_POLYNOMIAL_MODEL_HPP
