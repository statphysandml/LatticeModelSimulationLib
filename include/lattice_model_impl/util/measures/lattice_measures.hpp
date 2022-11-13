#ifndef LATTICEMODELIMPLEMENTATIONS_SYSTEM_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_SYSTEM_MEASURES_HPP

#include <mcmc/mcmc_simulation/measure_policy.hpp>
#include <mcmc/mcmc_simulation/util/random.hpp>
#include <param_helper/params.hpp>
#include <nlohmann/json.hpp>


namespace lm_impl {
    namespace util {
        namespace system_measures {
            template<typename SB>
            struct MeasureEnergyPolicy : public mcmc::measures::Measure<SB> {
            public:
                static auto compute_measure(const SB &system) {
                    return system.energy() / double(system.size());
                }

                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name() {
                    return "Energy";
                }
            };

            template<typename SB>
            struct MeasureDriftPolicy : public mcmc::measures::Measure<SB> {
            public:
                static auto compute_measure(const SB &system) {
                    return system.drift_term() / double(system.size());
                }

                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name() {
                    return "Drift";
                }
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP
