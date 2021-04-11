//
// 
//

#ifndef LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SUTWO_HPP
#define LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SUTWO_HPP


#include "../../mcmc_update_base.hpp"
#include "../../../site/site_models/onelink_su_two_model.hpp"

#include "../../../../external/Eigen/Dense"
#include "../../../../external/unsupported/Eigen/MatrixFunctions"

using namespace Eigen;
using namespace std::literals;


namespace lm_impl {
    namespace mcmc_update {

        template<typename ModelParameters>
        class LangevinUpdateSUTwo;


        template<typename ModelParameters>
        class LangevinUpdateSUTwoParameters : public MCMCUpdateBaseParameters {
        public:
            explicit LangevinUpdateSUTwoParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                             epsilon(get_entry<double>("epsilon", eps)),
                                                                             sqrt2epsilon(sqrt(2 * get_entry<double>(
                                                                                     "epsilon", eps))) {}

            explicit LangevinUpdateSUTwoParameters(
                    const double epsilon_
            ) : LangevinUpdateSUTwoParameters(json{
                    {"epsilon", epsilon_},
                    {"eps",     epsilon_}
            }) {}

            static std::string name() {
                return "LangevinUpdateSUTwo";
            }

            typedef LangevinUpdateSUTwo<ModelParameters> MCMCUpdate;

        private:
            friend class LangevinUpdateSUTwo<ModelParameters>;

            const double epsilon;
            const double sqrt2epsilon;
        };


        template<typename ModelParameters>
        class LangevinUpdateSUTwo
                : public MCMCUpdateBase<LangevinUpdateSUTwo<ModelParameters>, lm_impl::site_system::OneLinkSUTwoModelSampler> {
        public:
            explicit LangevinUpdateSUTwo(const LangevinUpdateSUTwoParameters<ModelParameters> &up_,
                                             typename ModelParameters::Model &model_)
                    : MCMCUpdateBase<LangevinUpdateSUTwo<ModelParameters>, lm_impl::site_system::OneLinkSUTwoModelSampler>(up_.eps), up(up_),
                      model(model_) {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            auto estimate_drift_term(const T site) {
                return model.get_drift_term(site);
            }

            template<typename T>
            T operator()(const T site) {

                std::complex<double> K1, K2, K3;
                std::tie(K1, K2, K3) = model.get_K_terms(site);
                std::complex<double> K_terms[3] = {K1, K2, K3};

                return update(site, &K_terms[3], up.epsilon, up.sqrt2epsilon);
            }

            template<typename T>
            T operator()(const T site, const double KMax, const double KExpectation) {

                std::complex<double> K1, K2, K3;
                std::tie(K1, K2, K3) = model.get_K_terms(site);
                std::complex<double> eps_drift_term[3] = {K1, K2, K3};
                double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
                return update(site, &eps_drift_term[3], epsilon, std::sqrt(2 * epsilon));
            }

        private:
            const LangevinUpdateSUTwoParameters<ModelParameters> &up;
            typename ModelParameters::Model &model;
            std::vector<double> epsilon;

            std::normal_distribution<double> normal;

            template<typename T>
            T update(T link, const std::complex<double> K_values[3], const double &epsilon, const double &sqrt2epsilon) {
                //T new_site(0);

                Matrix2cd pauli_1;
                pauli_1 << 0, 1.0,
                            1.0, 0;

                Matrix2cd pauli_2;
                pauli_2 << 0.0, -1i,
                            1i, 0.0;

                Matrix2cd pauli_3;
                pauli_3 << 1.0, 0.0,
                            0.0, -1.0;

                std::vector<Matrix2cd> pauli_matrices = {pauli_1, pauli_2, pauli_3};

                Matrix2cd M;
                M << 0.0, 0.0, 0.0, 0.0;
                for (int i = 0; i < 3; i++) {
                    M += 1i*pauli_matrices[i]*(epsilon*K_values[i] + sqrt2epsilon*normal(mcmc::util::gen));
                }

                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();

                Matrix2cd U;
                U << x0, x1, x2, x3;

                Matrix2cd R = M.exp();

                Matrix2cd U_new = R*U;

                T new_link(U_new(0,0), U_new(0,1), U_new(1,0), U_new(1,1));
                
                return new_link;
                
            }
        };
    }
}

#endif //LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SUTWO_HPP
