#ifndef MAIN_SITE_HPP
#define MAIN_SITE_HPP

#include <mcmc_simulation/header.hpp>
#include <param_helper/params.hpp>


#include "../util/measures/lattice_measures.hpp"


namespace lm_impl {
    namespace site_system {

        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        class SiteSystem
                : public mcmc::simulation::MeasureSystemBase<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>> {
        public:
            explicit SiteSystem(const json params) : mcmc::simulation::MeasureSystemBase<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>>(params) {
                sampler_ptr_ = std::make_unique<Sampler>(
                    mcmc::util::generate_parameter_class_json<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>, Sampler>(
                        *this, Sampler::name()));

                mcmc_model_ptr_ = std::make_unique<Model>(
                    mcmc::util::generate_parameter_class_json<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>, Model>(
                        *this, Model::name()));

                mcmc_method_ptr_ = std::make_unique<Method>(
                    mcmc::util::generate_parameter_class_json<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>, Method>(
                        *this, Method::name()));

                site_update_ptr_ = std::make_unique<UpdateDynamics>(
                    mcmc::util::generate_parameter_class_json<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>, UpdateDynamics>(
                        *this, UpdateDynamics::name()));

                // Needs to stay here since other site_update or update_formalism use site as reference
                initialize_site();

                sampler_ptr_->init(*this);
                mcmc_model_ptr_->init(*this);
                mcmc_method_ptr_->init(*this);
                site_update_ptr_->init(*this);
            }

            SiteSystem(
                    Sampler sampler=Sampler(),
                    Model mcmc_model=Model(),
                    Method mcmc_method=Method(),
                    UpdateDynamics lattice_update=UpdateDynamics(),
                    const std::vector<std::string> measures=std::vector<std::string> {}
            ) : SiteSystem(
                json {
                    {Sampler::name(), sampler.get_json()},
                    {Model::name(), mcmc_model.get_json()},
                    {Method::name(), mcmc_method.get_json()},
                    {UpdateDynamics::name(), lattice_update.get_json()},
                    {"measures", measures}
                }
            )
            {}

            typedef SiteSystem<T, Model, Method, UpdateDynamics, Sampler> System;

            void write_to_file(const std::string rel_config_path) override {
                std::string sampler_params_path = this->template get_entry<std::string>(
                    Sampler::name() + "_path", rel_config_path);
                sampler_ptr_->write_to_file(sampler_params_path);

                std::string mcmc_model_params_path = this->template get_entry<std::string>(
                    Model::name() + "_path", rel_config_path);
                mcmc_model_ptr_->write_to_file(mcmc_model_params_path);

                std::string mcmc_method_params_path = this->template get_entry<std::string>(
                    Method::name() + "_path", rel_config_path);
                mcmc_method_ptr_->write_to_file(mcmc_method_params_path);

                std::string site_update_params_path = this->template get_entry<std::string>(
                    UpdateDynamics::name() + "_path", rel_config_path);
                site_update_ptr_->write_to_file(site_update_params_path);

                json sampler_parameters = sampler_ptr_->get_json();
                this->template delete_entry(Sampler::name());

                json mcmc_model_parameters = mcmc_model_ptr_->get_json();
                this->template delete_entry(Model::name());

                json update_parameters = mcmc_method_ptr_->get_json();
                this->template delete_entry(Method::name());

                json site_update_parameters = site_update_ptr_->get_json();
                this->template delete_entry(UpdateDynamics::name());

                param_helper::params::Parameters::write_to_file(rel_config_path, this->name());

                this->template add_entry(Sampler::name(), sampler_parameters);
                this->template add_entry(Model::name(), mcmc_model_parameters);
                this->template add_entry(Method::name(), update_parameters);
                this->template add_entry(UpdateDynamics::name(), site_update_parameters);
            }

            param_helper::params::Parameters build_expanded_raw_parameters() const override {
                param_helper::params::Parameters parameters(this->params_);
                parameters.add_entry(Sampler::name(), sampler_ptr_->get_json());
                parameters.add_entry(Model::name(), mcmc_model_ptr_->get_json());
                parameters.add_entry(Method::name(), mcmc_method_ptr_->get_json());
                parameters.add_entry(UpdateDynamics::name(), site_update_ptr_->get_json());
                return parameters;
            }

            void initialize(std::string starting_mode) {
                // Needs to be called at the end so that generated objects (dynamics, mcmc_model, etc.) can already be used!
                this->generate_measures(this->measure_names());

                if (starting_mode == "hot")
                    site_ = sampler_ptr_->template random_state<T>();
                else
                    site_ = sampler_ptr_->template cold_state<T>();
            }

            void update_step(uint measure_interval) {
                site_update_ptr_->operator()(*this, measure_interval);
            }

            auto get_size() const {
                return 1;
            }

            auto at(int i) const {
                return site_;
            }

            auto &at(int i) {
                return site_;
            }

            auto get_system_representation() const {
                return site_;
            }

            auto &get_system_representation() {
                return site_;
            }

            void generate_measures(const std::vector<std::string> &measure_names) override {
                this->measurements_.clear();
                auto site_related_measures = generate_site_system_measures(this->measure_names());
                this->concat_measures(site_related_measures);

                auto sampler_related_measures = sampler_ptr_-> template generate_sampler_measures<
                    SiteSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(sampler_related_measures);

                auto mcmc_model_related_measures = mcmc_model_ptr_-> template generate_mcmc_model_measures<
                    SiteSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(mcmc_model_related_measures);

                auto mcmc_method_related_measures = mcmc_method_ptr_-> template generate_mcmc_method_measures<
                    SiteSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(mcmc_method_related_measures);

                auto site_update_related_measures = site_update_ptr_-> template generate_update_dynamics_measures<
                    SiteSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(site_update_related_measures);

                auto common_defined_measures = this->generate_systembase_measures(this->measure_names());
                this->concat_measures(common_defined_measures);
            }

            auto& get_sampler() const {
                return sampler_ptr_;
            }

            auto& get_sampler() {
                return sampler_ptr_;
            }

            auto& get_mcmc_model() const {
                return mcmc_model_ptr_;
            }

            auto& get_mcmc_model() {
                return mcmc_model_ptr_;
            }

            auto& get_mcmc_method() const {
                return mcmc_method_ptr_;
            }

            auto& get_mcmc_method() {
                return mcmc_method_ptr_;
            }

            auto& get_site_update() const {
                return site_update_ptr_;
            }

            auto& get_site_update() {
                return site_update_ptr_;
            }

            auto energy() const {
                return mcmc_model_ptr_->get_potential(site_);
            }

            auto drift_term() const {
                return mcmc_model_ptr_->get_drift_term(site_);
            }

            void normalize(T &site_elem) {
                site_elem = mcmc_model_ptr_->normalize_state(site_elem);
            }

            typedef T SiteType;

        protected:
            std::unique_ptr<Sampler> sampler_ptr_;
            std::unique_ptr<Model> mcmc_model_ptr_;
            std::unique_ptr<Method> mcmc_method_ptr_;
            std::unique_ptr<UpdateDynamics> site_update_ptr_;

            T site_;

            void initialize_site();

            std::vector<std::unique_ptr<mcmc::measures::Measure<SiteSystem>>>
            generate_site_system_measures(const std::vector<std::string> &measure_names);
        };


        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        void
        SiteSystem<T, Model, Method, UpdateDynamics, Sampler>::initialize_site() {
            site_ = sampler_ptr_->template random_state<T>();
        }

        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        std::ostream &operator<<(std::ostream &os,
                                 const SiteSystem<T, Model, Method, UpdateDynamics, Sampler> &site) {
            std::cout << site(0) << std::endl;
            return os;
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        std::vector<std::unique_ptr<mcmc::measures::Measure<SiteSystem<T, Model, Method, UpdateDynamics, Sampler>>>>
        SiteSystem<T, Model, Method, UpdateDynamics, Sampler>::generate_site_system_measures(
                const std::vector<std::string> &measure_names) {
            typedef SiteSystem<T, Model, Method, UpdateDynamics, Sampler> SiteSys;
            std::vector<std::unique_ptr<mcmc::measures::Measure<SiteSys>>> site_measures{};
            for (auto &measure_name :  measure_names)
                if (measure_name == "Energy")
                    site_measures.push_back(std::make_unique<util::system_measures::MeasureEnergyPolicy<SiteSys>>());
 //               else if (measure_name == "Drift")
 //                   site_measures.push_back(std::make_unique<util::system_measures::MeasureDriftPolicy<SiteSys>>());
            return site_measures;
        }

    }
}


#endif //MAIN_SITE_HPP
