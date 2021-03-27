//
// Created by lukas on 20.03.21.
//

#ifndef LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_ON_HPP
#define LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_ON_HPP


#include "../../mcmc_update_base.hpp"
#include "../../../lattice/lattice_models/on_model.hpp"


namespace lm_impl {
    namespace mcmc_update {

        template<typename ModelParameters>
        class LangevinUpdateON;


        template<typename ModelParameters>
        class LangevinUpdateONParameters : public MCMCUpdateBaseParameters {
        public:
            explicit LangevinUpdateONParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                             epsilon(get_entry<double>("epsilon", eps)),
                                                                             sqrt2epsilon(sqrt(2 * get_entry<double>(
                                                                                     "epsilon", eps))) {}

            explicit LangevinUpdateONParameters(
                    const double epsilon_
            ) : LangevinUpdateONParameters(json{
                    {"epsilon", epsilon_},
                    {"eps",     epsilon_}
            }) {}

            static std::string name() {
                return "LangevinUpdateON";
            }

            typedef LangevinUpdateON<ModelParameters> MCMCUpdate;

        private:
            friend class LangevinUpdateON<ModelParameters>;

            const double epsilon;
            const double sqrt2epsilon;
        };


        template<typename ModelParameters>
        class LangevinUpdateON
                : public MCMCUpdateBase<LangevinUpdateON<ModelParameters>, lm_impl::lattice_system::ONModelSampler> {
        public:
            explicit LangevinUpdateON(const LangevinUpdateONParameters<ModelParameters> &up_,
                                             typename ModelParameters::Model &model_)
                    : MCMCUpdateBase<LangevinUpdateON<ModelParameters>, lm_impl::lattice_system::ONModelSampler>(up_.eps), up(up_),
                      model(model_) {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            T estimate_drift_term(const T site) {
                return model.get_drift_term(site);
            }

            template<typename T>
            T estimate_drift_term(const T site, const std::vector<T *> neighbours) {
                return model.get_drift_term(site, neighbours);
            }

            template<typename T>
            T operator()(const T site) {
                return update(site, model.get_drift_term(site), up.epsilon, up.sqrt2epsilon);
            }

            template<typename T>
            T operator()(const T site, const std::vector<T *> neighbours) {
                return update(site, model.get_drift_term(site, neighbours), up.epsilon, up.sqrt2epsilon);
            }

            template<typename T>
            T operator()(const T site, const double KMax, const double KExpectation) {
                T eps_drift_term = model.get_drift_term(site);
                double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
                return update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
            }

            template<typename T>
            T
            operator()(const T site, const std::vector<T *> neighbours, const double KMax, const double KExpectation) {
                T eps_drift_term = model.get_drift_term(site, neighbours);
                double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
                return update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
            }

        private:
            const LangevinUpdateONParameters<ModelParameters> &up;
            typename ModelParameters::Model &model;
            std::vector<double> epsilon;

            std::normal_distribution<double> normal;

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++) {
                    new_site(i) = site(i) - epsilon * drift_term(i) + sqrt2epsilon * normal(mcmc::util::gen);
                }
                return model.normalize(new_site);
            }
        };
    }
}

#endif //LATTICEMODELSIMULATIONLIB_LANGEVIN_UPDATE_ON_HPP
