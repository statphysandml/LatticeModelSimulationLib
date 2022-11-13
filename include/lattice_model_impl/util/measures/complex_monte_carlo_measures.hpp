#ifndef COMPLEXMONTECARLO_COMPLEX_MONTE_CARLO_MEASURES_HPP
#define COMPLEXMONTECARLO_COMPLEX_MONTE_CARLO_MEASURES_HPP

#include <mcmc/mcmc_simulation/measure_policy.hpp>
#include <lattice_model_impl/site/site.hpp>
#include <lattice_model_impl/mcmc_method/dummy_mcmc_method.hpp>
#include <lattice_model_impl/update_dynamics/dummy_update_dynamics.hpp>


namespace lm_impl {
    namespace util {
        namespace complex_monte_carlo_measures {
            template<typename SB>
            struct MeasureDetailedBalanceAccuracyPolicy : public mcmc::measures::Measure<SB> {
            public:
                explicit MeasureDetailedBalanceAccuracyPolicy(const json &system_repr_params_,
                                                              const typename SB::SiteType &prev_site_) :
                        prev_site(prev_site_), system_repr(
                        site_system::SiteSystem<typename SB::SiteType, typename SB::ModelType, mcmc_method::DummyMCMCMethodParameters, update_dynamics::DummyUpdateDynamicsParameters>::from_json(
                                json{{"mcmc_method_params", {{"eps", 0.0}}},
                                     {"model_params_path",  system_repr_params_["model_params_path"]}})) {}

                std::string measure(const SB &system) override {
                    auto previous_site = prev_site;

                    system_repr() = previous_site;
                    auto previous_energy = system_repr.energy();
                    previous_site.real(system().real());
                    system_repr() = previous_site;
                    auto current_energy = system_repr.energy();

                    return std::to_string(current_energy.imag() - previous_energy.imag());
                }

                std::string name() {
                    return "DetailedBalanceAccuracy";
                }

                const typename SB::SiteType &prev_site;
                site_system::SiteSystem<typename SB::SiteType, typename SB::ModelType, mcmc_method::DummyMCMCMethodParameters, update_dynamics::DummyUpdateDynamicsParameters> system_repr;
            };

            template<typename SB>
            struct MeasureAbsoluteDetailedBalanceAccuracyPolicy : public mcmc::measures::Measure<SB> {
            public:
                explicit MeasureAbsoluteDetailedBalanceAccuracyPolicy(const json &system_repr_params_,
                                                                      const typename SB::SiteType &prev_site_) :
                        prev_site(prev_site_), system_repr(
                        site_system::SiteSystem<typename SB::SiteType, typename SB::ModelType, mcmc_method::DummyMCMCMethodParameters, update_dynamics::DummyUpdateDynamicsParameters>::from_json(
                                json{{"mcmc_method_params", {{"eps", 0.0}}},
                                     {"model_params_path",  system_repr_params_["model_params_path"]}})) {}

                std::string measure(const SB &system) override {
                    auto previous_site = prev_site;

                    system_repr() = previous_site;
                    auto previous_energy = system_repr.energy();
                    previous_site.real(system().real());
                    system_repr() = previous_site;
                    auto current_energy = system_repr.energy();

                    return std::to_string(std::abs(current_energy.imag() - previous_energy.imag()));
                }

                std::string name() {
                    return "AbsoluteDetailedBalanceAccuracy";
                }

                const typename SB::SiteType &prev_site;
                site_system::SiteSystem<typename SB::SiteType, typename SB::ModelType, mcmc_method::DummyMCMCMethodParameters, update_dynamics::DummyUpdateDynamicsParameters> system_repr;
            };

            template<typename SB>
            struct MeasureRealStepSizePolicy : public mcmc::measures::Measure<SB> {
            public:
                explicit MeasureRealStepSizePolicy(typename SB::SiteType &prev_site_) :
                        prev_site(prev_site_) {}

                std::string measure(const SB &system) override {
                    return std::to_string(std::abs(system().real() - prev_site.real()));
                }

                std::string name() {
                    return "RealStepSize";
                }

                typename SB::SiteType &prev_site;
            };

            template<typename SB>
            struct MeasureImagStepSizePolicy : public mcmc::measures::Measure<SB> {
            public:
                explicit MeasureImagStepSizePolicy(typename SB::SiteType &prev_site_) :
                        prev_site(prev_site_) {}

                std::string measure(const SB &system) override {
                    return std::to_string(std::abs(system().imag() - prev_site.imag()));
                }

                std::string name() {
                    return "ImagStepSize";
                }

                typename SB::SiteType &prev_site;
            };

            template<typename SB>
            struct MeasureComplexStepSizePolicy : public mcmc::measures::Measure<SB> {
            public:
                explicit MeasureComplexStepSizePolicy(typename SB::SiteType &prev_site_) :
                        prev_site(prev_site_) {}

                std::string measure(const SB &system) override {
                    return std::to_string(std::abs(system() - prev_site));
                }

                std::string name() {
                    return "ComplexStepSize";
                }

                typename SB::SiteType &prev_site;
            };
        }
    }
}


#endif //COMPLEXMONTECARLO_COMPLEX_MONTE_CARLO_MEASURES_HPP
