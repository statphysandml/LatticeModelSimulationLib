//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP


#include "mcmc_simulation/util/random.hpp"
#include "../mcmc_model_base.hpp"


namespace lm_impl {
    namespace lattice_system {
        class ONModel : public lm_impl::model::MCMCModelBase<ONModel> {
        public:
            explicit ONModel(const json params):
                MCMCModelBase(params),
                kappa_(set_kappa()),
                lambda_(set_lambda())
            {}

            ONModel(const double kappa_=1.0, const double lambda_=1.0):
                ONModel(json {
                    {"kappa", kappa_},
                    {"lambda", lambda_}
                })
            {}

            double set_kappa()
            {
                if(haskey("m0") and haskey("lambda0"))
                {
                    auto dim = get_entry<double>("dim");
                    double m0 = get_entry<double>("m0");
                    double lambda0 = get_entry<double>("lambda0") / 6.0;  // Adpatation to 1/4! instead of 1/4;
                    auto kappa = (-1.0 * (m0 + 2.0 * dim) + std::sqrt(std::pow(m0 + 2.0 * dim, 2.0) + 8 * lambda0)) / (4.0 * lambda0);
                    add_entry("kappa", kappa);
                }
                return get_entry<double>("kappa");
            }

            double set_lambda()
            {
                if(haskey("m0") and haskey("lambda0"))
                {
                    auto dim = get_entry<double>("dim");
                    double m0 = get_entry<double>("m0");
                    double lambda0 = get_entry<double>("lambda0")  / 6.0;  // Adpatation to 1/4! instead of 1/4;
                    auto kappa = (-1.0 * (m0 + 2.0 * dim) + std::sqrt(std::pow(m0 + 2.0 * dim, 2.0) + 8 * lambda0)) / (4.0 * lambda0);
                    add_entry("lambda", lambda0 * std::pow(kappa, 2.0));
                }
                return get_entry<double>("lambda");
            }

            static ONModel from_lambda_and_kappa(double lambda, double kappa)
            {
                return ONModel(json{
                        {"kappa", kappa},
                        {"lambda", lambda}
                });
            }

            static ONModel from_lambda_and_mass(double lambda0, double m0, double dim)
            {
                return ONModel(json{
                        {"dim", dim},
                        {"m0", m0},
                        {"lambda0", lambda0}
                });
            }

            // According to equation (3.16) from Smit QFT Lattice with the adaption the we choose 1/4! instead of 1/4
            template<typename T, typename T2=double_t>
            T2 get_potential(const T site, const std::vector<T*> neighbours) const {
                double potential = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    potential += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                potential = 2.0 * kappa_ * potential - site_sq - lambda_ * pow(site_sq - 1.0, 2.0);
                return -1.0 * potential;
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours) const {
                double potential = 0;
                for(size_t i = 0; i < neighbours.size(); i += 2) {
                    potential += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                potential = 2.0 * kappa_ * potential - site_sq - lambda_ * pow(site_sq - 1.0, 2.0);
                return -1.0 * potential;
            }

            template<typename T>
            T get_drift_term(const T site, const std::vector<T*> neighbours) const {
                auto drift = T(0);
                for(size_t i = 0; i < neighbours.size(); i++) {
                    drift += *neighbours[i];
                }
                double site_sq = site * site;
                drift = -2.0 * kappa_ * drift + 2.0 * site + 4.0 * lambda_ * site * (site_sq - 1.0);
                return drift;
            }

        private:
            double kappa_;
            double lambda_;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP
