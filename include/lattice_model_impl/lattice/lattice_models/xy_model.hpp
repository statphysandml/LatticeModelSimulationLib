//
// Created by lukas on 31.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP


#include "../mcmc_model_base.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"


namespace lm_impl {
    namespace lattice_system {

        template<typename SB>
        struct MeasureXYMagnetization : public mcmc::measures::Measure<SB> {
        public:
            static auto compute_measure(const SB &system) {
                auto sum_cos = decltype(system[0]){0};
                auto sum_sin = decltype(system[0]){0};

                for (uint i = 0; i < system.size(); i++) {
                    sum_cos += std::cos(system[i]);
                    sum_sin += std::sin(system[i]);
                }
                return std::pow(sum_cos / system.size(), 2) + std::pow(sum_sin / system.size(), 2);
            }

            std::string measure(const SB &system) override {
                return std::to_string(compute_measure(system));
            }

            std::string name() {
                return "XYMagnetization";
            }
        };

        class XYModel : public lm_impl::model::MCMCModelBase<XYModel> {
        public:
            explicit XYModel(const json params):
                MCMCModelBase(params),
                beta_(get_entry<double>("beta", 0.5)),
                J_(get_entry<double>("J", 1.0)),
                h_(get_entry<double>("h", 0.0))
            {}

            explicit XYModel(double beta=0.5, double J=1.0, double h=0.0) : XYModel(json{
                    {"beta", beta},
                    {"J",    J},
                    {"h",    h}
            }) {}

            template<typename T>
            T normalize(T state) {
                state = std::fmod(state, 2 * M_PI);
                if (state < 0) {
                    state += 2 * M_PI;
                }
                return state;
            }

            template<typename T>
            T get_potential(const T site, const std::vector<T*> neighbours) const
            {
                double S = 0;
                for (size_t i = 0; i < neighbours.size(); i ++) {
                    S += J_ * std::cos(site - *neighbours[i]);
                }
                S += h_ * std::cos(site);
                return -1.0 * beta_ * S;
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours) const
            {
                double S = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S += J_ * std::cos(site - *neighbours[i]);
                }
                S += h_ * std::cos(site);
                return -1.0 * beta_ * S;
            }

            template<typename T>
            T get_drift_term(const T site, const std::vector<T*> neighbours) const
            {
                double S = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S += J_ * std::sin(site - *neighbours[i]) +
                         J_ * std::sin(site - *neighbours[i + 1]);
                }
                S += h_ * std::sin(site);
                return beta_ * S;
            }

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>
            generate_model_measures(const SB &system) {
                auto measure_names = system.measure_names();

                std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>> measurements{};
                for (auto &measure_name :  measure_names)
                    if (measure_name == "XYMagnetization")
                        measurements.push_back(std::make_unique<MeasureXYMagnetization <SB>>());
                return measurements;
            }

        private:
            double beta_;
            double J_;
            double h_;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_XY_MODEL_HPP
