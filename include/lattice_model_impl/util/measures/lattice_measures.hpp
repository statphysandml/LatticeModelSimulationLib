//
// Created by lukas on 10.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP

#include <mcmc_simulation/measure_policy.hpp>
#include <mcmc_simulation/util/random.hpp>
#include <param_helper/params.hpp>
#include <param_helper/json.hpp>


namespace lm_impl {
    namespace util {
        namespace lattice_system_model_measures {
            template<typename SB>
            struct MeasureEnergyPolicy : public mcmc::measures::Measure<SB> {
            public:
                std::string measure(const SB &system) override {
                    return std::to_string(system.energy() / double(system.size()));
                }

                std::string name() {
                    return "Energy";
                }
            };

            template<typename SB>
            struct MeasureDriftPolicy : public mcmc::measures::Measure<SB> {
            public:
                std::string measure(const SB &system) override {
                    return std::to_string(system.drift_term() / double(system.size()));
                }

                std::string name() {
                    return "Drift";
                }
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_LATTICE_MEASURES_HPP
