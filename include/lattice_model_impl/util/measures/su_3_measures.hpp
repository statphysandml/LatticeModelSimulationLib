//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SU3_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_SU3_MEASURES_HPP

#include <mcmc_simulation/measure_policy.hpp>
#include <mcmc_simulation/util/random.hpp>
#include <param_helper/params.hpp>
#include <param_helper/json.hpp>


namespace lm_impl {
    namespace util {
        namespace su_3_measures {

            template<typename SB>
            struct MeasurePolyakovLoopPolicy : public mcmc::measures::Measure<SB> {
            public:
                std::string measure(const SB &system) override {
                    auto model = system.get_model();
                    auto site = system.get_system_representation();
                    auto P = model.get_P(site);

                    return std::to_string(P);
                }

                std::string name() {
                    return "Polyakov_Loop";
                }
            };

            template<typename SB>
            struct MeasureConjugatePolyakovLoopPolicy : public mcmc::measures::Measure<SB> {
            public:
                std::string measure(const SB &system) override {
                    auto model = system.get_model();
                    auto site = system.get_system_representation();
                    auto P_inv = model.get_P_inv(site);

                    return std::to_string(P_inv);
                }

                std::string name() {
                    return "Conjugate_Polyakov_Loop";
                }
            };

            template<typename SB>
            struct MeasurePolyakovLoop_realPolicy : public mcmc::measures::Measure<SB> {
            public:
                std::string measure(const SB &system) override {
                    auto model = system.get_model();
                    auto site = system.get_system_representation();
                    auto P_real = std::real(model.get_P(site));

                    return std::to_string(P_real);
                }

                std::string name() {
                    return "Polyakov_Loop_real";
                }
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SU3_MEASURES_HPP
