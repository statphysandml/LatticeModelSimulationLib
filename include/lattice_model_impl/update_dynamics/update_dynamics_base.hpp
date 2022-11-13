#ifndef LATTICEMODELIMPLEMENTATIONS_UPDATE_DYNAMICS_BASE_HPP
#define LATTICEMODELIMPLEMENTATIONS_UPDATE_DYNAMICS_BASE_HPP


#include <param_helper/params.hpp>
#include <mcmc/mcmc_simulation/measure_policy.hpp>

namespace lm_impl {
    namespace update_dynamics {
        template<typename UpdateFormalism, typename... Args>
        auto update_system_site(UpdateFormalism &f, Args &&... args) {
            return f(std::forward<Args>(args)...);
        }

        template<typename UpdateFormalism, typename... Args>
        void global_system_update(UpdateFormalism &f, Args &&... args) {
            f(std::forward<Args>(args)...);
        }

        template<typename UpdateDynamics>
        class UpdateDynamicsBase : public param_helper::params::Parameters{
        public:
            using Parameters::Parameters;

            void write_to_file(const std::string rel_root_dir) {
                Parameters::write_to_file(rel_root_dir, UpdateDynamicsBase::name());
            }

            static std::string name() {
                return "update_dynamics";
            }

            template<typename System>
            void init(System &system) {
                return update_dynamics().initialize(system);
            }

            template<typename System>
            void operator()(System &system, uint measure_interval = 1) {
                return update_dynamics().update(system, measure_interval);
            }

            template<typename SB>
            std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>
            generate_update_dynamics_measures(const SB &system) {
                // auto measure_names = system_parameters.measure_names();
                return std::vector<std::unique_ptr<mcmc::measures::Measure<SB>>>{};
            }

        protected:
            template<typename System>
            void initialize(System &system) {}

        private:
            UpdateDynamics &update_dynamics() {
                return *static_cast<UpdateDynamics*>(this);
            }

            const UpdateDynamics &update_dynamics() const {
                return *static_cast<const UpdateDynamics*>(this);
            }
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_UPDATE_DYNAMICS_BASE_HPP
