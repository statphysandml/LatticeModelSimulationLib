//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_LINK_LATTICE_MODEL_HPP
#define MAIN_LINK_LATTICE_MODEL_HPP

#include "../lattice/mcmc_model_base.hpp"

#include "../util/measures/link_lattice_model_measures.hpp"

namespace lm_impl {
    namespace link_lattice_system {

        template <typename Model>
        class LinkLatticeModel : public lm_impl::model::MCMCModelBase<Model>
        {
        public:
            using lm_impl::model::MCMCModelBase<Model>::MCMCModelBase;

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>
            generate_model_measures(const SB &system) {
                auto measure_names = system.measure_names();
                std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>> link_lattice_measures{};
                for (auto &measure_name :  measure_names)
                    if (measure_name == "PolyakovLoop")
                        link_lattice_measures.push_back(std::make_unique<util::link_lattice_system_model_measures::MeasurePolyakovLoopPolicy<SB>>(
                                system.get_dimensions(), system.get_elems_per_site()));
                    else if (measure_name == "AveragePlaquetteAction")
                        link_lattice_measures.push_back(std::make_unique<util::link_lattice_system_model_measures::MeasureAveragePlaquetteActionPolicy<SB>>(
                                system.get_dimension()));
                return link_lattice_measures;
            }

        protected:
            template<typename T>
            T calc_A (const std::vector<T*> neighbours, bool both_orientations=true) const {
                T A("null");
                for(uint i = 0; i < neighbours.size(); i += 6)
                {
                    A += (*neighbours[i]) * T(*neighbours[i + 1]).adjungate() * T(*neighbours[i + 2]).adjungate();
                    if(both_orientations)
                        A += T(*neighbours[i + 3]).adjungate() * T(*neighbours[ i + 4]).adjungate() * (*neighbours[i + 5]);
                }
                return A;
            }
        };

    }
}


#endif //MAIN_LINK_LATTICE_MODEL_HPP
