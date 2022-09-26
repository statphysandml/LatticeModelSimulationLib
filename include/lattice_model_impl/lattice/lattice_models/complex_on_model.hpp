//
// Created by lukas on 12.01.21.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_ON_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_ON_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"
#include "../mcmc_model_base.hpp"


namespace lm_impl {
    namespace lattice_system {
        class ComplexONModel : public lm_impl::model::MCMCModelBase<ComplexONModel> {
        public:
            explicit ComplexONModel(const json params) :
                MCMCModelBase(params),
                kappa_(std::complex<double> {get_entry<double>("kappa_real", 1.0),
                                             get_entry<double>("kappa_imag", 0.0)}),
                lambda_(std::complex<double> {get_entry<double>("lambda_real", 1.0),
                                              get_entry<double>("lambda_imag", 0.0)})
            {}

            explicit ComplexONModel(double kappa_real=1.0, double kappa_imag=0.0,
                double lambda_real=1.0, double lambda_imag=0.0) : ComplexONModel(json{
                    {"kappa_real", kappa_real},
                    {"kappa_imag", kappa_imag},
                    {"lambda_real", lambda_real},
                    {"lambda_imag", lambda_imag}
            })
            {}

            // According to equation (3.16) from Smit QFT Lattice
            template<typename T, typename T2=std::complex<double_t>>
            T2 get_potential(const T site, const std::vector<T*> neighbours) const {
                std::complex<double> potential = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    potential += site * (*neighbours[i]);
                }
                std::complex<double> site_sq = site * site;
                potential = 2.0 * kappa_ * potential - site_sq - lambda_ * pow(site_sq - 1.0, 2.0);
                return -1.0 * potential;
            }

            template<typename T, typename T2=std::complex<double_t>>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours) const {
                std::complex<double> potential = 0;
                for(size_t i = 0; i < neighbours.size(); i += 2) {
                    potential += site * (*neighbours[i]);
                }
                std::complex<double> site_sq = site * site;
                potential = 2.0 * kappa_ * potential - site_sq - lambda_ * pow(site_sq - 1.0, 2.0);
                return -1.0 * potential;
            }

            template<typename T>
            T get_drift_term(const T site, const std::vector<T*> neighbours) const {
                T drift_term(0);
                for(size_t i = 0; i < neighbours.size(); i++) {
                    drift_term += (*neighbours[i]);
                }
                drift_term = 2.0 * kappa_ * drift_term + (-2.0) * site + (-4.0) * site * lambda_ * (site * site - 1.0);
                return -1.0 * drift_term;
            }

        private:
            std::complex<double> kappa_;
            std::complex<double> lambda_;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_ON_MODEL_HPP
