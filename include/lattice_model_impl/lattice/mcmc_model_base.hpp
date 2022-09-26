#ifndef MAIN_MODEL_HPP
#define MAIN_MODEL_HPP


#include <param_helper/params.hpp>
#include <mcmc_simulation/measure_policy.hpp>


namespace lm_impl {
    namespace model {
        template<typename MCMCModel>
        class MCMCModelBase : public param_helper::params::Parameters {
        public:
            using Parameters::Parameters;

            void write_to_file(const std::string rel_root_dir) {
                Parameters::write_to_file(rel_root_dir, "mcmc_model_params");
            }

            static const std::string name() {
                return "mcmc_model";
            }

            template<typename System>
            void init(System &site) {
                return lattice_model().initialize(site);
            }

            template<typename T>
            T normalize_state(T state) {
                return lattice_model().normalize(state);
            }

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>
            generate_mcmc_model_measures(const SB &system) {
                // auto measure_names = system_parameters.get_measures();
                return std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>{};
            }

        protected:
            template<typename System>
            void initialize(System &system) {}

            template<typename T>
            T normalize(T state) {
                return state;
            }

        private:
            MCMCModel &lattice_model() {
                return *static_cast<MCMCModel*>(this);
            }

            const MCMCModel &system_update() const {
                return *static_cast<const MCMCModel*>(this);
            }
        };

    }
}

#endif //MAIN_MODEL_HPP
