#ifndef LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SU3_HPP
#define LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SU3_HPP

#include "langevin_update_base.hpp"

//#include "../../mcmc_update_base.hpp"
#include "../../mcmc_method_base.hpp"
#include "../../../site/site_models/onelink_su_3_model.hpp"

#include "../../../../external/Eigen/Dense"
#include "../../../../external/unsupported/Eigen/MatrixFunctions"

using namespace Eigen;
using namespace std::literals;


namespace lm_impl {
    namespace mcmc_method {
        template<typename Model>
        class ComplexLangevinSU3ModelUpdate : public LangevinUpdateBase<ComplexLangevinSU3ModelUpdate<Model>, Model> {
        public:
            using LangevinUpdateBase<ComplexLangevinSU3ModelUpdate<Model>, Model>::LangevinUpdateBase;

            template<typename T>
            T update(T link, const std::array<std::complex<double>, 8> K_values, const double &epsilon, const double &sqrt2epsilon) {

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
                    M += 1i*lambdas[i]*(epsilon*K_values[i] + sqrt2epsilon*this->normal_(mcmc::util::g_gen));
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

                T new_link(U_new(0,0), U_new(0,1), U_new(0,2),
                            U_new(1, 0), U_new(1,1), U_new(1,2),
                            U_new(2,0), U_new(2, 1), U_new(2, 2));
                
                return new_link;
                
            }
        };
    }
}

#endif //LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_SU3_HPP
