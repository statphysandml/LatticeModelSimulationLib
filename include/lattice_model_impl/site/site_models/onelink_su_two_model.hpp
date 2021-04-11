//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_OneLink_SU_TWO_MODEL_HPP
#define MAIN_OneLink_SU_TWO_MODEL_HPP

#include "../site_model.hpp"
//#include "../link_lattice_model.hpp"
#include <math.h>

#include "../../../external/Eigen/Dense"
#include "../../../external/unsupported/Eigen/MatrixFunctions"

using namespace Eigen;
using namespace std::literals;


namespace lm_impl {
    namespace site_system {


        class OneLinkSUTwoModel;

        class OneLinkSUTwoModelParameters : public SiteModelParameters {
        public:
            explicit OneLinkSUTwoModelParameters(const json params_) : SiteModelParameters(params_),
                                                                beta(get_entry<std::complex<double>>("beta")),
                                                                mu(get_entry<std::complex<double>>("mu")),
                                                                kappa(get_entry<std::complex<double>>("kappa"))
                                                                {}

            explicit OneLinkSUTwoModelParameters(std::complex<double> beta_,
                                          std::complex<double> mu_, 
                                          std::complex<double> kappa_) : OneLinkSUTwoModelParameters(
            json{
                    {"beta", beta_},
                    {"mu", mu_},
                    {"kappa", kappa_}
            }) {}

            const static std::string name() {
                return "SUTwoModel";
            }

            static uint N() {
                return 2;
            };

            typedef OneLinkSUTwoModel Model;

        private:
            friend class OneLinkSUTwoModel;

            const std::complex<double> beta;
            const std::complex<double> mu;
            const std::complex<double> kappa;
        };


        class OneLinkSUTwoModel : public SiteModel<OneLinkSUTwoModel> {
        public:
            explicit OneLinkSUTwoModel(const OneLinkSUTwoModelParameters &mp_) : mp(mp_) {}

            template<typename T>
            std::complex<double> get_potential(T link) {
                
                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();

                Matrix2cd U;
                U << x0, x1, x2, x3;

                return -mp.beta / (2.0*OneLinkSUTwoModelParameters::N()) * (U.trace() + U.inverse().trace());
            }

            template<typename T>
            std::complex<double> get_P(T link) {
                
                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();

                Matrix2cd U;
                U << x0, x1, x2, x3;

                return 1.0 / OneLinkSUTwoModelParameters::N() * U.trace();
            }

            template<typename T>
            std::complex<double> get_P_inv(T link) {
                
                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();

                Matrix2cd U;
                U << x0, x1, x2, x3;

                return 1.0 / OneLinkSUTwoModelParameters::N() * U.inverse().trace();
            }

            template<typename T>
            std::complex<double> get_n(T link) {
                
                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();

                Matrix2cd U;
                U << x0, x1, x2, x3;

                Matrix2cd U_inv = U.inverse();

                auto P = 1.0/OneLinkSUTwoModelParameters::N() * U.trace();
                auto P_inv = 1.0/OneLinkSUTwoModelParameters::N() * U_inv.trace();

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
            std::tuple<T2, T2, T2> get_K_terms(T link) {

                Matrix2cd pauli_1;
                pauli_1 << 0, 1.0,
                            1.0, 0;

                Matrix2cd pauli_2;
                pauli_2 << 0.0, -1i,
                            1i, 0.0;

                Matrix2cd pauli_3;
                pauli_3 << 1.0, 0.0,
                            0.0, -1.0;

                Matrix2cd I;
                pauli_3 << 1.0, 0.0,
                            0.0, 1.0;

                auto x0 = link.x0();
                auto x1 = link.x1();
                auto x2 = link.x2();
                auto x3 = link.x3();

                Matrix2cd U;
                U << x0, x1, x2, x3;
                Matrix2cd U_inv = U.inverse();

                auto P = 1/3.0 * U.trace();
                auto P_inv = 1/3.0 * U_inv.trace();

                std::complex<double> mu = mp.mu;
                std::complex<double> kappa = mp.kappa;

                auto M_q = 1.0+3.0*kappa*exp(mu)*P + 3.0*pow(kappa, 2.0)*exp(2.0*mu)*P_inv + pow(kappa, 3.0)*exp(3.0*mu);
                auto M_qbar = 1.0+3.0*kappa*exp(-mu)*P_inv + 3.0*pow(kappa, -2.0)*exp(2.0*mu)*P + pow(kappa, 3.0)*exp(-3.0*mu);

                auto D1P = 1i/3.0 * (pauli_1*U).trace();
                auto D2P = 1i/3.0 * (pauli_2*U).trace();
                auto D3P = 1i/3.0 * (pauli_3*U).trace();

                auto D1P_inv = -1i/3.0 * (pauli_1*U_inv).trace();
                auto D2P_inv = -1i/3.0 * (pauli_2*U_inv).trace();
                auto D3P_inv = -1i/3.0 * (pauli_3*U_inv).trace();

                auto K1B = mp.beta/2.0 *(D1P + D1P_inv);
                auto K1F = 3.0/M_q * (kappa*exp(mu)*D1P + pow(kappa, 2.0)*exp(2.0*mu)*D1P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D1P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D1P);

                auto K2B = mp.beta/2.0 *(D2P + D2P_inv);
                auto K2F = 3.0/M_q * (kappa*exp(mu)*D2P + pow(kappa, 2.0)*exp(2.0*mu)*D2P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D2P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D2P);

                auto K3B = mp.beta/2.0 *(D3P + D3P_inv);
                auto K3F = 3.0/M_q * (kappa*exp(mu)*D3P + pow(kappa, 2.0)*exp(2.0*mu)*D3P_inv)
                            + 3.0/M_qbar * (kappa*exp(-mu)*D3P_inv + pow(kappa, 2.0)*exp(-2.0*mu)*D3P);
                
                return std::make_tuple(K1B + K1F, K2B + K2F, K3B + K3F);
            }


        private:
            const OneLinkSUTwoModelParameters &mp;
        };

        struct OneLinkSUTwoModelSampler
        {
            OneLinkSUTwoModelSampler(const double eps_) : eps(eps_)
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
                return "SUTwoModelSampler";
            }

            const double eps;
        };
    }
}

#endif //MAIN_OneLink_SU_TWO_MODEL_HPP
