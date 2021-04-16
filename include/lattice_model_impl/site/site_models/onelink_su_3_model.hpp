//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_OneLink_SU_3_MODEL_HPP
#define MAIN_OneLink_SU_3_MODEL_HPP

#include "../site_model.hpp"
#include <math.h>

#include "../../../external/Eigen/Dense"
#include "../../../external/unsupported/Eigen/MatrixFunctions"

using namespace Eigen;
using namespace std::literals;


namespace lm_impl {
    namespace site_system {



        template<typename SB>
        struct MeasureOneLinkPlaquette : public mcmc::common_measures::MeasurePolicy<SB> {
        public:
            std::string measure(const SB &system) override {
                auto model = system.get_model();
                auto site = system.get_system_representation();
                auto P = model.get_P(site);

                return std::to_string(P);
            }

            std::string name() {
                return "OneLinkPlaquette";
            }
        };

        template<typename SB>
        struct MeasureOneLinkInversePlaquette : public mcmc::common_measures::MeasurePolicy<SB> {
        public:
            std::string measure(const SB &system) override {
                auto model = system.get_model();
                auto site = system.get_system_representation();
                auto P_inv = model.get_P_inv(site);

                return std::to_string(P_inv);
            }

            std::string name() {
                return "OneLinkInversePlaquette";
            }
        };


        class OneLinkSU3Model;

        class OneLinkSU3ModelParameters : public SiteModelParameters {
        public:
            explicit OneLinkSU3ModelParameters(const json params_) : SiteModelParameters(params_),

                                                                beta(std::complex<double> {get_entry<double>("beta_real", 0.0),
                                                                                    get_entry<double>("beta_imag", 0.0)}),
                                                                mu(std::complex<double> {get_entry<double>("mu_real", 0.0),
                                                                                    get_entry<double>("mu_imag", 0.0)}),
                                                                kappa(std::complex<double> {get_entry<double>("kappa_real", 0.0),
                                                                                    get_entry<double>("kappa_imag", 0.0)})
            {}

            explicit OneLinkSU3ModelParameters(double beta_real_, double beta_imag_, double mu_real_,
                                                      double mu_imag_, double kappa_real_, double kappa_imag_)
                    : OneLinkSU3ModelParameters(json{
                    {"beta_real", beta_real_},
                    {"beta_imag", beta_imag_},
                    {"mu_real", mu_real_},
                    {"mu_imag", mu_imag_},
                    {"kappa_real", kappa_real_},
                    {"kappa_imag", kappa_imag_},
            }) {}

            const static std::string name() {
                return "OneLinkSU3Model";
            }

            static uint N() {
                return 3;
            };

            typedef OneLinkSU3Model Model;

        private:
            friend class OneLinkSU3Model;

            const std::complex<double> beta;
            const std::complex<double> mu;
            const std::complex<double> kappa;
        };


        class OneLinkSU3Model : public SiteModel<OneLinkSU3Model> {
        public:
            explicit OneLinkSU3Model(const OneLinkSU3ModelParameters &mp_) : mp(mp_) {}

            template<typename T>
            std::complex<double> get_potential(T link) {
                
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

                return -mp.beta / (2.0*OneLinkSU3ModelParameters::N()) * (U.trace() + U.inverse().trace());
            }

            template<typename T>
            std::complex<double> get_P(T link) {
                
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

                return 1.0 / OneLinkSU3ModelParameters::N() * U.trace();
            }

            template<typename T>
            std::complex<double> get_P_inv(T link) {
                
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

                return 1.0 / OneLinkSU3ModelParameters::N() * U.inverse().trace();
            }

            template<typename T>
            std::complex<double> get_n(T link) {
                
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

                Matrix3cd U_inv = U.inverse();

                std::complex<double> P = 1.0/OneLinkSU3ModelParameters::N() * U.trace();
                std::complex<double> P_inv = 1.0/OneLinkSU3ModelParameters::N() * U_inv.trace();

                std::complex<double> mu = mp.mu;
                std::complex<double> kappa = mp.kappa;

                auto M_q = 1.0+3.0*kappa*exp(mu)*P + 3.0*pow(kappa, 2.0)*exp(2.0*mu)*P_inv + pow(kappa, 3.0)*exp(3.0*mu);
                auto M_qbar = 1.0+3.0*kappa*exp(-mu)*P_inv + 3.0*pow(kappa, -2.0)*exp(2.0*mu)*P + pow(kappa, 3.0)*exp(-3.0*mu);

                auto n = 3.0/M_q * (kappa*exp(mu)*P + 2.0*pow(kappa, 2.0)*exp(2.0*mu)*P_inv + pow(kappa, 3.0)*exp(3.0*mu))
                            + 3.0/M_qbar * (kappa*exp(-mu)*P_inv + 2.0*pow(kappa, 2.0)*exp(-2.0*mu)*P + pow(kappa, 3.0)*exp(-3.0*mu));

                return n;
            }




            template<typename T>
            std::complex<double> get_drift_term(T link) {
                return 0.0;
            }


            template<typename T, typename T2=std::complex<double>>
            std::array<std::complex<double>, 8> get_K_terms(T link) {

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
                Matrix3cd U_inv = U.inverse();

                auto P = 1/3.0 * U.trace();
                auto P_inv = 1/3.0 * U_inv.trace();

                std::complex<double> mu = mp.mu;
                std::complex<double> kappa = mp.kappa;

                auto M_q = 1.0+3.0*kappa*exp(mu)*P + 
                            3.0*pow(kappa, 2.0)*exp(2.0*mu)*P_inv + pow(kappa, 3.0)*exp(3.0*mu);
                
                auto M_qbar = 1.0+3.0*kappa*exp(-mu)*P_inv + 
                            3.0*pow(kappa, 2.0)*exp(-2.0*mu)*P + pow(kappa, 3.0)*exp(-3.0*mu);

                auto D1P = 1i/3.0 * (lambda_1*U).trace();
                auto D2P = 1i/3.0 * (lambda_2*U).trace();
                auto D3P = 1i/3.0 * (lambda_3*U).trace();
                auto D4P = 1i/3.0 * (lambda_4*U).trace();
                auto D5P = 1i/3.0 * (lambda_5*U).trace();
                auto D6P = 1i/3.0 * (lambda_6*U).trace();
                auto D7P = 1i/3.0 * (lambda_7*U).trace();
                auto D8P = 1i/3.0 * (lambda_8*U).trace();

                auto D1P_inv = -1i/3.0 * (U_inv*lambda_1).trace();
                auto D2P_inv = -1i/3.0 * (U_inv*lambda_2).trace();
                auto D3P_inv = -1i/3.0 * (U_inv*lambda_3).trace();
                auto D4P_inv = -1i/3.0 * (U_inv*lambda_4).trace();
                auto D5P_inv = -1i/3.0 * (U_inv*lambda_5).trace();
                auto D6P_inv = -1i/3.0 * (U_inv*lambda_6).trace();
                auto D7P_inv = -1i/3.0 * (U_inv*lambda_7).trace();
                auto D8P_inv = -1i/3.0 * (U_inv*lambda_8).trace();

                auto K1B = mp.beta/2.0 *(D1P + D1P_inv);
                auto K1F = 3.0/M_q * (kappa*exp(mu)*D1P + pow(kappa, 2.0)*exp(2.0*mu)*D1P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D1P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D1P);

                auto K2B = mp.beta/2.0 *(D2P + D2P_inv);
                auto K2F = 3.0/M_q * (kappa*exp(mu)*D2P + pow(kappa, 2.0)*exp(2.0*mu)*D2P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D2P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D2P);

                auto K3B = mp.beta/2.0 *(D3P + D3P_inv);
                auto K3F = 3.0/M_q * (kappa*exp(mu)*D3P + pow(kappa, 2.0)*exp(2.0*mu)*D3P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D3P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D3P);

                auto K4B = mp.beta/2.0 *(D4P + D4P_inv);
                auto K4F = 3.0/M_q * (kappa*exp(mu)*D4P + pow(kappa, 2.0)*exp(2.0*mu)*D4P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D4P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D4P);

                auto K5B = mp.beta/2.0 *(D5P + D5P_inv);
                auto K5F = 3.0/M_q * (kappa*exp(mu)*D5P + pow(kappa, 2.0)*exp(2.0*mu)*D5P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D5P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D5P);

                auto K6B = mp.beta/2.0 *(D6P + D6P_inv);
                auto K6F = 3.0/M_q * (kappa*exp(mu)*D6P + pow(kappa, 2.0)*exp(2.0*mu)*D6P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D6P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D6P);

                auto K7B = mp.beta/2.0 *(D7P + D7P_inv);
                auto K7F = 3.0/M_q * (kappa*exp(mu)*D7P + pow(kappa, 2.0)*exp(2.0*mu)*D7P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D7P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D7P);

                auto K8B = mp.beta/2.0 *(D8P + D8P_inv);
                auto K8F = 3.0/M_q * (kappa*exp(mu)*D8P + pow(kappa, 2.0)*exp(2.0*mu)*D8P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D8P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D8P);

                std::array<std::complex<double>, 8> values = { K1B + K1F, K2B + K2F, K3B + K3F, K4B + K4F,
                                                                K5B + K5F, K6B + K6F, K7B + K7F, K8B + K8F };   
            return values;
            }

        template<typename SB, typename SBP>
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>>
            generate_model_measures(const SBP &system_parameters) {
                auto measure_names = system_parameters.get_measures();

                std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SB>>> measures{};
                for (auto &measure_name :  measure_names) {
                    if (measure_name == "OneLinkPlaquette")
                        measures.push_back(std::make_unique<MeasureOneLinkPlaquette <SB>>());
                    if (measure_name == "OneLinkInversePlaquette")
                        measures.push_back(std::make_unique<MeasureOneLinkInversePlaquette <SB>>());
                };
                return measures;
            }


        private:
            const OneLinkSU3ModelParameters &mp;
        };

        struct OneLinkSU3ModelSampler
        {
            OneLinkSU3ModelSampler(const double eps_) : eps(eps_)
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
                return "OneLinkSU3ModelSampler";
            }

            const double eps;
        };
    }
}

#endif //MAIN_OneLink_SU_TWO_MODEL_HPP
