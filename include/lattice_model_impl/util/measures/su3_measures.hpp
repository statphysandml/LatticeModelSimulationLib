#ifndef LATTICEMODELIMPLEMENTATIONS_SU3_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_SU3_MEASURES_HPP

#include <mcmc/mcmc_simulation/measure_policy.hpp>
#include <mcmc/mcmc_simulation/util/random.hpp>
#include <param_helper/params.hpp>
#include <nlohmann/json.hpp>


namespace lm_impl {
    namespace util {
        namespace su_3_measures {

            template<typename SB>
            struct MeasurePolyakovLoopPolicy : public mcmc::measures::Measure<SB> {
            public:
                static auto compute_measure(const SB &system)
                {
                    auto model_ptr_ = system.get_mcmc_model().get();
                    auto site = system.get_system_representation();
                    auto P = model_ptr_->get_P(site);

                    return P;
                }
 
                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name() override {
                        return "Polyakov_Loop";
                }
            };

            template<typename SB>
            struct MeasureConjugatePolyakovLoopPolicy : public mcmc::measures::Measure<SB> {
            public:
                static auto compute_measure(const SB &system) {
                    auto model_ptr_ = system.get_mcmc_model().get();
                    auto site = system.get_system_representation();
                    auto P_inv = model_ptr_->get_P_inv(site);

                    return P_inv;
                }

                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name() override {
                    return "Conjugate_Polyakov_Loop";
                }
            };

            template<typename SB>
            struct MeasurePolyakovLoop_realPolicy : public mcmc::measures::Measure<SB> {
            public:
                static auto compute_measure(const SB &system) {
                    auto model_ptr_ = system.get_mcmc_model().get();
                    auto site = system.get_system_representation();
                    auto P_real = std::real(model_ptr_->get_P(site));

                    return P_real;
                }

                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name() override {
                    return "Polyakov_Loop_real";
                }
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SU3_MEASURES_HPP
