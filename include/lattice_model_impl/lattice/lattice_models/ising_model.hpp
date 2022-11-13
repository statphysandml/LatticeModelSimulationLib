#ifndef MAIN_ISING_MODEL_HPP
#define MAIN_ISING_MODEL_HPP

#include <mcmc/mcmc_simulation/util/random.hpp>
#include <nlohmann/json.hpp>

#include <lattice_model_impl/lattice/mcmc_model_base.hpp>
#include <lattice_model_impl/sampler/sampler_base.hpp>


namespace lm_impl {
    namespace lattice_system {
        class IsingModel : public lm_impl::model::MCMCModelBase<IsingModel>
        {
        public:
            explicit IsingModel(const json params):
                MCMCModelBase(params),
                beta_(get_entry<double>("beta", 0.4)),
                J_(get_entry<double>("J", 1.0)),
                h_(get_entry<double>("h", 0.0))
            {}

            IsingModel(double beta=0.4, double J=1.0, double h=0.0) : IsingModel(json{
                    {"beta", beta},
                    {"J", J},
                    {"h", h}
            })
            {}

            static const std::string type() {
                return "IsingModel";
            }

            template<typename T, typename T2=double_t>
            T2 get_potential(const T site, const std::vector<T*> neighbours) const
            {
                double coupling = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    coupling += *neighbours[i];
                }
                return  -1.0 * beta_ * site * (J_ * coupling + h_);
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours) const
            {
                double coupling = 0;
                // Only neighbours in positive direction
                for(size_t i = 0; i < neighbours.size(); i += 2) {
                    coupling += *neighbours[i];
                }
                return  -1.0 * beta_ * site * (J_ * coupling + h_);
            }

        private:
            double beta_;
            double J_;
            double h_;
        };

        struct IsingModelSampler : public lm_impl::sampler::SamplerBase<IsingModelSampler>
        {
            explicit IsingModelSampler(json params):
                SamplerBase<IsingModelSampler>(params)
            {
                uniint_ = std::uniform_int_distribution<int>(0, 1);
            }

            explicit IsingModelSampler() : IsingModelSampler(json{})
            {}

            static const std::string type() {
                return "IsingModelSampler";
            }

            template<typename T>
            T cold_sample() {
                return 1;
            }

            template<typename T>
            T random_sample()
            {
                return 2 * uniint_(mcmc::util::random::g_gen) - 1;
            }

            template<typename T>
            T propose_sample(T site)
            {
                return -1 * site;
            }

            double get_eps() const
            {
                return 0.0;
            }

            std::uniform_int_distribution<int> uniint_;
        };
    }
}



#endif //MAIN_ISING_MODEL_HPP
