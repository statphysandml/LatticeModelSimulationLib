#ifndef LATTICEMODELIMPLEMENTATIONS_MCMC_METHOD_BASE_HPP
#define LATTICEMODELIMPLEMENTATIONS_MCMC_METHOD_BASE_HPP


#include <param_helper/params.hpp>
#include <mcmc_simulation/util/random.hpp>
#include "../sampler/dummy_sampler.hpp"


namespace lm_impl {
    namespace mcmc_method {

        template<typename MCMCMethod>
        class MCMCMethodBase : public param_helper::params::Parameters {
        public:
            explicit MCMCMethodBase(const json params) :
                Parameters(params)
            {}

            MCMCMethodBase() : MCMCMethodBase(json{})
            {}

            void write_to_file(const std::string rel_root_dir) {
                Parameters::write_to_file(rel_root_dir, "mcmc_method_params");
            }

            static std::string name() {
                return "mcmc_method";
            }

            template<typename System>
            void init(System &site) {
                return mcmc_method().initialize(site);
            }

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>
            generate_mcmc_method_measures(const SB &system) {
                // auto measure_names = system_parameters.get_measures();
                return std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>{};
            }

        protected:
            template<typename System>
            void initialize(System &system) {}

        private:
            MCMCMethod &mcmc_method() {
                return *static_cast<MCMCMethod*>(this);
            }

            const MCMCMethod &mcmc_method() const {
                return *static_cast<const MCMCMethod*>(this);
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_MCMC_METHOD_BASE_HPP
