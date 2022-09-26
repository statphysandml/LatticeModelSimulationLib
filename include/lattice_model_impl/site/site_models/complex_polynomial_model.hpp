#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP


#include "../../lattice/mcmc_model_base.hpp"
#include "mcmc_simulation/util/random.hpp"


namespace lm_impl {
    namespace site_system {
        class ComplexPolynomialModel : public lm_impl::model::MCMCModelBase<ComplexPolynomialModel> {
        public:
            explicit ComplexPolynomialModel(const json params):
                MCMCModelBase(params),
                lambda_(std::complex<double> {get_entry<double>("lambda_real", 1.0), get_entry<double>("lambda_imag", 0.0)}),
                sigma_(std::complex<double> {get_entry<double>("sigma_real", 1.0), get_entry<double>("sigma_imag", 1.0)}),
                h_(std::complex<double> {get_entry<double>("h_real", 0.0), get_entry<double>("h_imag", 0.0)})
            {}

            explicit ComplexPolynomialModel(double lambda_real=1.0, double lambda_imag=0.0, double sigma_real=1.0,
                double sigma_imag=1.0, double h_real=0.0, double h_imag=0.0):
                ComplexPolynomialModel(json{
                    {"lambda_real", lambda_real},
                    {"lambda_imag", lambda_imag},
                    {"sigma_real", sigma_real},
                    {"sigma_imag", sigma_imag},
                    {"h_real", h_real},
                    {"h_imag", h_imag},
            }) {}

            std::complex<double> get_potential(const std::complex<double> site) const {
                return 0.5 * sigma_ * std::pow(site, 2) + 0.25 * lambda_ * std::pow(site, 4) + h_ * site;
            }

            std::complex<double> get_drift_term(const std::complex<double> site) const {
                return sigma_ * site + lambda_ * std::pow(site, 3) + h_;
            }

        private:
            std::complex<double> lambda_;
            std::complex<double> sigma_;
            std::complex<double> h_;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_POLYNOMIAL_MODEL_HPP
