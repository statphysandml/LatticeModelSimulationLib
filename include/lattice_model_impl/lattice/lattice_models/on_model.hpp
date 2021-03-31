//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP

#include "mcmc_simulation/util/random.hpp"
#include "../lattice_model.hpp"

namespace lm_impl {
    namespace lattice_system {

        class ONModel;

        class ONModelParameters : public LatticeModelParameters {
        public:
            explicit ONModelParameters(const json params_) : LatticeModelParameters(params_),
                                                                kappa(set_kappa()),
                                                                lambda(set_lambda())
            {}

            double set_kappa()
            {
                if(haskey("m0") and haskey("lambda0"))
                {
                    auto dim = get_entry<double>("dim");
                    double m0_ = get_entry<double>("m0");
                    double lambda0_ = get_entry<double>("lambda0") / 6.0;  // Adpatation to 1/4! instead of 1/4;
                    auto kappa_ = (-1.0 * (m0_ + 2.0 * dim) + std::sqrt(std::pow(m0_ + 2.0 * dim, 2.0) + 8 * lambda0_)) / (4.0 * lambda0_);
                    add_entry("kappa", kappa_);
                }
                return get_entry<double>("kappa");
            }

            double set_lambda()
            {
                if(haskey("m0") and haskey("lambda0"))
                {
                    auto dim = get_entry<double>("dim");
                    double m0_ = get_entry<double>("m0");
                    double lambda0_ = get_entry<double>("lambda0")  / 6.0;  // Adpatation to 1/4! instead of 1/4;
                    auto kappa_ = (-1.0 * (m0_ + 2.0 * dim) + std::sqrt(std::pow(m0_ + 2.0 * dim, 2.0) + 8 * lambda0_)) / (4.0 * lambda0_);
                    add_entry("lambda", lambda0_ * std::pow(kappa_, 2.0));
                }
                return get_entry<double>("lambda");
            }

            static ONModelParameters from_lambda_and_kappa(double lambda_, double kappa_)
            {
                return ONModelParameters(json{
                        {"kappa", kappa_},
                        {"lambda", lambda_}
                });
            }

            static ONModelParameters from_lambda_and_mass(double lambda0_, double m0_, double dim_)
            {
                return ONModelParameters(json{
                        {"dim", dim_},
                        {"m0", m0_},
                        {"lambda0", lambda0_}
                });
            }

            const static std::string name() {
                return "ONModel";
            }

            typedef ONModel Model;

        private:
            friend class ONModel;

            const double kappa;
            const double lambda;
        };

        class ONModel : public LatticeModel<ONModel> {
        public:
            explicit ONModel(const ONModelParameters &mp_) : mp(mp_) {}

            // According to equation (3.16) from Smit QFT Lattice with the adaption the we choose 1/4! instead of 1/4
            template<typename T, typename T2=double_t>
            T2 get_potential(const T site, const std::vector<T*> neighbours)
            {
                double potential = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    potential += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                potential = 2.0 * mp.kappa * potential - site_sq - mp.lambda * pow(site_sq - 1.0, 2.0);
                return -1.0 * potential;
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours)
            {
                double potential = 0;
                for(size_t i = 0; i < neighbours.size(); i += 2) {
                    potential += site * (*neighbours[i]);
                }
                double site_sq = site * site;
                potential = 2.0 * mp.kappa * potential - site_sq - mp.lambda * pow(site_sq - 1.0, 2.0);
                return -1.0 * potential;
            }

            template<typename T>
            T get_drift_term(const T site, const std::vector<T *> neighbours) {
                auto drift = T(0);
                for(size_t i = 0; i < neighbours.size(); i++) {
                    drift += *neighbours[i];
                }
                double site_sq = site * site;
                drift = -2.0 * mp.kappa * drift + 2.0 * site + 4.0 * mp.lambda * site * (site_sq - 1.0);
                return drift;
            }

        private:
            const ONModelParameters &mp;
        };

        // Overload GaussianSampler instead?

        struct ONModelSampler
        {
            ONModelSampler(const double eps_) : eps(eps_)
            {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            T random_state() {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++)
                    new_site(i) += std::sqrt(2 * eps) * normal(mcmc::util::gen);
                return new_site;
            }

            template<typename T>
            T propose_state(T site) {
                return site + random_state<T>();
            }

            double get_eps() const
            {
                return eps;
            }

            const static std::string name() {
                return "ONModelSampler";
            }

            const double eps;
            std::normal_distribution<double> normal;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_ON_MODEL_HPP
