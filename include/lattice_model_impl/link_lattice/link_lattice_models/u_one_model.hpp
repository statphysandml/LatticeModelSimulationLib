//
// Created by lukas on 04.11.19.
//

#ifndef MAIN_U_ONE_MODEL_HPP
#define MAIN_U_ONE_MODEL_HPP

#include "../link_lattice_model.hpp"

namespace lm_impl {
    namespace link_lattice_system {


        class UOneModel;

        class UOneModelParameters : public LinkLatticeModelParameters {
        public:
            explicit UOneModelParameters(const json params_) : LinkLatticeModelParameters(params_),
                                                               beta(get_entry<double>("beta"))
            {}

            explicit UOneModelParameters(double beta_) : UOneModelParameters(json{
                    {"beta", beta_}
            }) {}

            const static std::string name() {
                return "UOneModel";
            }

            static uint N() {
                return 1;
            };

            typedef UOneModel Model;

        private:
            friend class UOneModel;

            const double beta;
        };


        class UOneModel : public LinkLatticeModel<UOneModel> {
        public:
            explicit UOneModel(const UOneModelParameters &mp_) : mp(mp_) {}

            template<typename T, typename T2=double_t>
            T2 get_potential(const T link, const std::vector<T *> neighbours) {
                T A = calc_A(neighbours);
                return mp.beta / UOneModelParameters::N() * (neighbours.size() / 3.0 - UOneModelParameters::N() *
                                           std::real((link * A).trace()));
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T link, const std::vector<T*> neighbours) {
                T A = calc_A(neighbours, false);
                return mp.beta / UOneModelParameters::N() * (neighbours.size() / (2.0 * 3.0) - UOneModelParameters::N() *
                                           std::real((link * A).trace()));
            }

            template<typename T, typename T2=double_t>
            T2 get_drift_term(const T link, const std::vector<T *> neighbours) {
                return T2(0.0);
            }

        private:
            const UOneModelParameters &mp;
        };

        struct UOneModelSampler
        {
            UOneModelSampler(const double eps_) : eps(eps_)
            {}

            template<typename T>
            T random_state() {
                return T("random");
            }

            template<typename T>
            T propose_state(T link) {
                return link * T(eps);
            }

            double get_eps() const
            {
                return eps;
            }

            const static std::string name() {
                return "UOneModelSampler";
            }

            const double eps;
        };
    }
}

#endif //MAIN_U_ONE_MODEL_HPP
