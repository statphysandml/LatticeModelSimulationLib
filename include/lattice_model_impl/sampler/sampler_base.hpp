#ifndef LATTICEMODELIMPLEMENTATIONS_SAMPLER_BASE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SAMPLER_BASE_HPP


#include <param_helper/params.hpp>
#include <mcmc/mcmc_simulation/measure_policy.hpp>
#include <mcmc/mcmc_simulation/util/random.hpp>

#include <lattice_model_impl/sampler/dummy_sampler.hpp>


namespace lm_impl {
    namespace sampler {

        template<typename Sampler>
        class SamplerBase : public param_helper::params::Parameters {
        public:
            explicit SamplerBase(const json params) :
                Parameters(params)
            {}

            SamplerBase() : SamplerBase(json{})
            {}

            void write_to_file(const std::string rel_root_dir) {
                Parameters::write_to_file(rel_root_dir, SamplerBase::name());
            }

            static std::string name() {
                return "sampler";
            }

            template<typename System>
            void init(System &site) {
                return sampler().initialize(site);
            }

            template<typename T>
            T random_state() {
                return sampler().template random_sample<T>();
            };

            template<typename T>
            T cold_state() {
                return sampler().template cold_sample<T>();
            };

            template<typename T>
            T propose_state(T site)
            {
                return sampler().template propose_sample<T>(site);
            }

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>
            generate_sampler_measures(const SB &system) {
                // auto measure_names = system_parameters.get_measures();
                return std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>{};
            }

        protected:
            template<typename System>
            void initialize(System &system) {}

        private:
            Sampler &sampler() {
                return *static_cast<Sampler*>(this);
            }

            const Sampler &sampler() const {
                return *static_cast<const Sampler*>(this);
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_SAMPLER_BASE_HPP
