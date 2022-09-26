//
// Created by lukas on 12.01.21.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_LINK_LATTICE_MODEL_MEASURES_HPP
#define LATTICEMODELIMPLEMENTATIONS_LINK_LATTICE_MODEL_MEASURES_HPP

#include "mcmc_simulation/measure_policy.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/params.hpp"
#include "param_helper/json.hpp"


namespace lm_impl {
    namespace util {
        namespace link_lattice_system_model_measures {
            template<typename SB>
            struct MeasurePolyakovLoopPolicy: public mcmc::measures::Measure< SB > {
            public:
                MeasurePolyakovLoopPolicy(const std::vector<int> dimensions, const uint elem_per_site) : dimensions_(dimensions), elem_per_site_(elem_per_site)
                {}

                auto compute_measure(const SB &system) {
                    // Hasn't been tested!!
                    std::complex<double> polyakov_loop = 0;
                    for(uint i = 0; i < system.size(); i += elem_per_site_ * dimensions_[0])
                    {
                        typename SB::SiteType group_elem = system[i];
                        for(uint tau = elem_per_site_; tau < elem_per_site_ * dimensions_[0]; tau += elem_per_site_)
                            group_elem = group_elem * system[i + tau];
                        polyakov_loop += group_elem.trace();
                    }
                    return polyakov_loop;
                }

                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name()
                {
                    return "PolyakovLoop";
                }

                const std::vector<int> dimensions_; // Different dimensions
                const uint elem_per_site_; // Number of elements per site
            };

            template<typename SB>
            struct MeasureAveragePlaquetteActionPolicy : public mcmc::measures::Measure<SB> {
            public:
                MeasureAveragePlaquetteActionPolicy(const int dimension) :
                    n_plaquettes_per_link_(dimension - 1) // Only in positive direction, since energy also computed only in positive direction
                {}

                // Needs to be devided by the inverse temperature in a postprocessing step
                auto compute_measure(const SB &system) {
                    return system.energy() / n_plaquettes_per_link_  / double(system.size());
                }

                std::string measure(const SB &system) override {
                    return std::to_string(compute_measure(system));
                }

                std::string name() {
                    return "AveragePlaquetteAction";
                }

                const int n_plaquettes_per_link_;
            };
        }
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_LINK_LATTICE_MODEL_MEASURES_HPP
