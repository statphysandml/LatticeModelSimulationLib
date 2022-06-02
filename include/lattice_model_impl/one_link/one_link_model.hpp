//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_ONE_LINK_MODEL_HPP
#define MAIN_ONE_LINK_MODEL_HPP

#include "../lattice/mcmc_model_base.hpp"

#include "../util/measures/su_3_measures.hpp"

namespace lm_impl {
    namespace one_link_system {

        template <typename Model>
        class OneLinkModel : public lm_impl::model::MCMCModelBase<Model>
        {
        public:
            using lm_impl::model::MCMCModelBase<Model>::MCMCModelBase;

            /*template<typename SB>
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
            }*/


            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>> 
            generate_model_measures(const SB &system) {
                auto measure_names = system.measure_names();

                std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>> su_3_measures{}; //measures{};
                for (auto &measure_name :  measure_names) {
                    if (measure_name == "PolyakovLoop")
                        su_3_measures.push_back(std::make_unique<util::su_3_measures::MeasurePolyakovLoopPolicy<SB>>());
                    if (measure_name == "ConjugatePolyakovLoop")
                        su_3_measures.push_back(std::make_unique<util::su_3_measures::MeasureConjugatePolyakovLoopPolicy<SB>>());
                        //measures.push_back(std::make_unique<MeasureOneLinkInversePlaquette <SB>>());
                    if (measure_name == "PolyakovLoop_real")
                        su_3_measures.push_back(std::make_unique<util::su_3_measures::MeasurePolyakovLoop_realPolicy<SB>>());
                        //measures.push_back(std::make_unique<MeasureOneLinkPlaquette_real <SB>>());
                };
                return su_3_measures;
            }

        protected:
            
        };

    }
}


#endif //MAIN_ONE_LINK_MODEL_HPP
