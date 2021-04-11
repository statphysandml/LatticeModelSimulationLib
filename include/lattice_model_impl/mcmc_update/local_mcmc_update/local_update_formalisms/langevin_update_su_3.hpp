//
// 
//

#ifndef LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SU3_HPP
#define LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SU3_HPP


#include "../../mcmc_update_base.hpp"
#include "../../../site/site_models/onelink_su_3_model.hpp"

#include "/home/adrian/bachelor_project/anharmonicoscillators/Hosak/MyProject/include/Eigen/Dense"
#include "/home/adrian/bachelor_project/anharmonicoscillators/Hosak/MyProject/include/unsupported/Eigen/MatrixFunctions"

using namespace Eigen;
using namespace std::literals;


namespace lm_impl {
    namespace mcmc_update {

        template<typename ModelParameters>
        class LangevinUpdateSU3;


        template<typename ModelParameters>
        class LangevinUpdateSU3Parameters : public MCMCUpdateBaseParameters {
        public:
            explicit LangevinUpdateSU3Parameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                             epsilon(get_entry<double>("epsilon", eps)),
                                                                             sqrt2epsilon(sqrt(2 * get_entry<double>(
                                                                                     "epsilon", eps))) {}

            explicit LangevinUpdateSU3Parameters(
                    const double epsilon_
            ) : LangevinUpdateSU3Parameters(json{
                    {"epsilon", epsilon_},
                    {"eps",     epsilon_}
            }) {}

            static std::string name() {
                return "LangevinUpdateSU3";
            }

            typedef LangevinUpdateSU3<ModelParameters> MCMCUpdate;

        private:
            friend class LangevinUpdateSU3<ModelParameters>;

            const double epsilon;
            const double sqrt2epsilon;
        };


        template<typename ModelParameters>
        class LangevinUpdateSU3
                : public MCMCUpdateBase<LangevinUpdateSU3<ModelParameters>, lm_impl::site_system::OneLinkSU3ModelSampler> {
        public:
            explicit LangevinUpdateSU3(const LangevinUpdateSU3Parameters<ModelParameters> &up_,
                                             typename ModelParameters::Model &model_)
                    : MCMCUpdateBase<LangevinUpdateSU3<ModelParameters>, lm_impl::site_system::OneLinkSU3ModelSampler>(up_.eps), up(up_),
                      model(model_) {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            auto estimate_drift_term(const T site) {
                return model.get_drift_term(site);
            }

            template<typename T>
            T operator()(const T site) {

                std::complex<double> K1, K2, K3, K4, K5, K6, K7, K8;
                std::tie(K1, K2, K3, K4, K5, K6, K7, K8) = model.get_K_terms(site);
                std::complex<double> K_terms[8] = {K1, K2, K3, K4, K5, K6, K7, K8};

                return update(site, &K_terms[8], up.epsilon, up.sqrt2epsilon);
            }

            template<typename T>
            T operator()(const T site, const double KMax, const double KExpectation) {

                std::complex<double> K1, K2, K3, K4, K5, K6, K7, K8;
                std::tie(K1, K2, K3, K4, K5, K6, K7, K8) = model.get_K_terms(site);
                std::complex<double> eps_drift_term[8] = {K1, K2, K3, K4, K5, K6, K7, K8};
                double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
                return update(site, &eps_drift_term[8], epsilon, std::sqrt(2 * epsilon));
            }

        private:
            const LangevinUpdateSU3Parameters<ModelParameters> &up;
            typename ModelParameters::Model &model;
            std::vector<double> epsilon;

            std::normal_distribution<double> normal;

            template<typename T>
            T update(T link, const std::complex<double> K_values[3], const double &epsilon, const double &sqrt2epsilon) {
                //T new_site(0);

                Matrix3cd lambda_1;
                lambda_1 << 0, 1.0, 0,
                            1.0, 0, 0,
                            0, 0, 0;

                Matrix3cd lambda_2;
                lambda_2 << 0, -1i, 0,
                            1i, 0, 0,
                            0, 0, 0;

                Matrix3cd lambda_3;
                lambda_3 <<  1.0, 0, 0,
                            0, -1.0, 0,
                            0, 0, 0;

                Matrix3cd lambda_4;
                lambda_4 <<  0, 0, 1.0,
                            0, 0, 0,
                            1.0, 0, 0;

                Matrix3cd lambda_5;
                lambda_5 <<  0, 0, -1i,
                            0, 0, 0,
                            1i, 0, 0;

                Matrix3cd lambda_6;
                lambda_6 <<  0, 0, 0,
                            0, 0, 1.0,
                            0, 1.0, 0;

                Matrix3cd lambda_7;
                lambda_7 <<  0, 0, 0,
                            0, 0, -1i,
                            0, 1i, 0;

                Matrix3cd lambda_8;
                lambda_8 <<  1.0, 0, 0,
                            0, 1.0, 0,
                            0, 0, -2.0;
                lambda_8 = 1/sqrt(3) * lambda_8;

                std::vector<Matrix3cd> lambdas = {lambda_1, lambda_2, lambda_3, lambda_4,
                                                    lambda_5, lambda_6, lambda_7, lambda_8};

                Matrix3cd M;
                M << 0, 0, 0, 0, 0, 0, 0, 0, 0;
                for (int i = 0; i < 8; i++) {
                    M += 1i*lambdas[i]*(epsilon*K_values[i] + sqrt2epsilon*normal(mcmc::util::gen));
                }

                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();
                auto x4 = link.x4();
                auto x5 = link.x5();
                auto x6 = link.x6();
                auto x7 = link.x7();
                auto x8 = link.x8();

                Matrix3cd U;
                U << x0, x1, x2, x3, x4, x5, x6, x7, x8;

                Matrix3cd R = M.exp();

                Matrix3cd U_new = R*U;

                T new_link(U_new(0,0), U_new(0,1), U_new(0,2), U_new(1, 0), U_new(1,1), U_new(1,2), U_new(2,0), U_new(2, 1), U_new(2, 2));
                
                return new_link;
                
            }
        };
    }
}

#endif //LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SU3_HPP
