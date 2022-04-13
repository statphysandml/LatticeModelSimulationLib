#ifndef MAIN_LATTICE_HPP
#define MAIN_LATTICE_HPP


#include <mcmc_simulation/header.hpp>
#include <param_helper/params.hpp>


#include "../util/measures/lattice_measures.hpp"


namespace lm_impl {
    namespace lattice_system {

        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        class LatticeSystem : public mcmc::simulation::MeasureSystemBase<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>> {
        public:
            explicit LatticeSystem(const json params) : mcmc::simulation::MeasureSystemBase<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>>(params)
            {
                dimensions_ = this->template get_entry<std::vector<int>>("dimensions", std::vector<int>{4, 4});
                lattice_action_type_ = this->template get_entry<std::string>("lattice_action_type", "nearest_neighbour");

                n_sites_ = 1;
                dim_mul_.push_back(n_sites_);
                for (auto dim: dimensions_) {
                    n_sites_ *= dim;
                    dim_mul_.push_back(n_sites_);
                }
                dimension_ = dimensions_.size();

                if(lattice_action_type_ == "nearest_neighbour")
                    elem_per_site_ = 1;
                else
                    elem_per_site_ = dimension_;

                sampler_ = Sampler(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>, Sampler>(
                        *this, Sampler::name()));

                mcmc_model_ = Model(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>, Model>(
                        *this, Model::name()));

                mcmc_method_ = Method(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>, Method>(
                        *this, Method::name()));

                lattice_update_ = UpdateDynamics(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>, UpdateDynamics>(
                        *this, UpdateDynamics::name()));

                initialize_lattice();
                if (lattice_action_type_ == "plaquette")
                    set_plaquette_neighbours();
                else
                    set_nearest_neighbours();

                sampler_.init(*this);
                mcmc_model_.init(*this);
                mcmc_method_.init(*this);
                lattice_update_.init(*this);
            }
            
            LatticeSystem(
                    Sampler sampler=Sampler(),
                    Model mcmc_model=Model(),
                    Method mcmc_method=Method(),
                    UpdateDynamics lattice_update=UpdateDynamics(),
                    const std::vector<int> dimensions=std::vector<int> {4, 4},
                    const std::string lattice_action_type="nearest_neighbour",
                    const std::vector<std::string> measures=std::vector<std::string> {}
            ) : LatticeSystem(
                json {
                    {Sampler::name(), sampler.get_json()},
                    {Model::name(), mcmc_model.get_json()},
                    {Method::name(), mcmc_method.get_json()},
                    {UpdateDynamics::name(), lattice_update.get_json()},
                    {"dimensions", dimensions},
                    {"lattice_action_type", lattice_action_type},
                    {"measures", measures}
                }
            )
            {}

            typedef LatticeSystem<T, Model, Method, UpdateDynamics, Sampler> System;

            void write_to_file(const std::string rel_config_path) override {
                std::string sampler_params_path = this->template get_entry<std::string>(
                    Sampler::name() + "_path", rel_config_path);
                sampler_.write_to_file(sampler_params_path);

                std::string mcmc_model_params_path = this->template get_entry<std::string>(
                    Model::name() + "_path", rel_config_path);
                mcmc_model_.write_to_file(mcmc_model_params_path);

                std::string mcmc_method_params_path = this->template get_entry<std::string>(
                    Method::name() + "_path", rel_config_path);
                mcmc_method_.write_to_file(mcmc_method_params_path);

                std::string lattice_update_params_path = this->template get_entry<std::string>(
                    UpdateDynamics::name() + "_path", rel_config_path);
                lattice_update_.write_to_file(lattice_update_params_path);

                json sampler_parameters = sampler_.get_json();
                this->template delete_entry(Sampler::name());

                json mcmc_model_parameters = mcmc_model_.get_json();
                this->template delete_entry(Model::name());

                json update_parameters = mcmc_method_.get_json();
                this->template delete_entry(Method::name());

                json lattice_update_parameters = lattice_update_.get_json();
                this->template delete_entry(UpdateDynamics::name());

                param_helper::params::Parameters::write_to_file(rel_config_path, this->name());

                this->template add_entry(Sampler::name(), sampler_parameters);
                this->template add_entry(Model::name(), mcmc_model_parameters);
                this->template add_entry(Method::name(), update_parameters);
                this->template add_entry(UpdateDynamics::name(), lattice_update_parameters);
            }

            param_helper::params::Parameters build_expanded_raw_parameters() const override {
                param_helper::params::Parameters parameters(this->params_);
                parameters.add_entry(Sampler::name(), sampler_.get_json());
                parameters.add_entry(Model::name(), mcmc_model_.get_json());
                parameters.add_entry(Method::name(), mcmc_method_.get_json());
                parameters.add_entry(UpdateDynamics::name(), lattice_update_.get_json());
                return parameters;
            }

            const std::vector<int> &get_dimensions() const {
                return dimensions_;
            }

            const auto &get_dimension() const {
                return dimension_;
            }

            const uint &get_elems_per_site() const {
                return elem_per_site_;
            }

            void initialize(std::string starting_mode) {
                // Needs to be called at the end so that generated objects (dynamics, mcmc_model, etc.) can already be used!
                this->generate_measures(this->measure_names());

                if(starting_mode == "hot")
                    for (auto &site : lattice_)
                        site = sampler_.template random_state<T>();
                else // starting_mode == "cold"
                    for (auto &site : lattice_)
                        site = sampler_.template cold_state<T>();
            }

            void update_step(uint measure_interval) {
                lattice_update_(*this, measure_interval);
            }

            // Returns the total number of elements of the lattice - (equals the total number of sites if elem_per_site=1)
            auto get_size() const {
                return n_sites_ * elem_per_site_;
            }

            auto at(int i) const {
                return lattice_[i];
            }

            auto &at(int i) {
                return lattice_[i];
            }

            auto get_system_representation() const {
                return lattice_;
            }

            auto &get_system_representation() {
                return lattice_;
            }

            auto get_elem_per_site() const {
                return elem_per_site_;
            }

            auto &neighbours_at(int i) {
                return neighbours_[i];
            }

            auto neighbours_at(int i) const {
                return neighbours_[i];
            }

            auto &get_neighbours() {
                return neighbours_;
            }

            auto get_neighbours() const {
                return neighbours_;
            }

            void generate_measures(const std::vector<std::string> &measure_names) override {
                this->measurements_.clear();
                auto lattice_related_measures = generate_lattice_system_measures(this->measure_names());
                this->concat_measures(lattice_related_measures);

                auto sampler_related_measures = sampler_. template generate_sampler_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(sampler_related_measures);

                auto mcmc_model_related_measures = mcmc_model_. template generate_mcmc_model_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(mcmc_model_related_measures);

                auto mcmc_method_related_measures = mcmc_method_. template generate_mcmc_method_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(mcmc_method_related_measures);

                auto lattice_update_related_measures = lattice_update_. template generate_update_dynamics_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>>(*this);
                this->concat_measures(lattice_update_related_measures);

                auto common_defined_measures = this->generate_systembase_measures(this->measure_names());
                this->concat_measures(common_defined_measures);
            }

            auto& get_sampler() const {
                return sampler_;
            }

            auto& get_sampler() {
                return sampler_;
            }

            auto& get_mcmc_model() const {
                return mcmc_model_;
            }

            auto& get_mcmc_model() {
                return mcmc_model_;
            }

            auto& get_mcmc_method() const {
                return mcmc_method_;
            }

            auto& get_mcmc_method() {
                return mcmc_method_;
            }

            auto& get_lattice_update() const {
                return lattice_update_;
            }

            auto& get_lattice_update() {
                return lattice_update_;
            }

            auto energy() const {
                decltype(mcmc_model_.get_energy_per_lattice_elem(lattice_[0], neighbours_[0])) energy(0);
                for (uint i = 0; i < get_size(); i++) {
                    energy += mcmc_model_.get_energy_per_lattice_elem(lattice_[i], neighbours_[i]);
                }
                return energy;
            }

            auto drift_term() const {
                std::cerr << "Drift term computation needs to be implemented here" << std::endl;
                std::exit(EXIT_FAILURE);
                decltype(mcmc_model_.get_potential(lattice_[0], neighbours_[0])) drift_term(0);
                for (uint i = 0; i < get_size(); i++) {
                    // ToDo: How can this be integrated?
                    // drift_term += mcmc_model.get_drift_term(lattice_[i], neighbours[i]);
                }
                return drift_term;
            }

            void normalize(std::vector<T> &lattice_grid) {
                for (auto &elem : lattice_grid)
                    elem = mcmc_model_.normalize_state(elem);
            }

            typedef T SiteType;

        protected:
            Sampler sampler_;
            Model mcmc_model_;
            Method mcmc_method_;
            UpdateDynamics lattice_update_;

            std::string lattice_action_type_;

            uint n_sites_; // Total number of sites
            uint elem_per_site_; // Number of elements per site (is equal to dimension), but only for the link mcmc_model
            std::vector<int> dimensions_; // Different dimensions
            std::vector<int> dim_mul_; // Accumulated different dimensions (by product)
            int dimension_; // Number of dimensions
            
            std::vector<T> lattice_;
            std::vector<std::vector<T*> > neighbours_;

            void initialize_lattice();

            int neigh_dir(int n, int d, bool dir, int mu) const;

            void set_nearest_neighbours();

            void set_plaquette_neighbours();

            std::vector<std::unique_ptr<mcmc::measures::Measure<LatticeSystem>>>
            generate_lattice_system_measures(const std::vector<std::string> &measure_names);
        };


        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        void LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>::initialize_lattice() {
            lattice_ = std::vector<T>(get_size(), T(0));
        }


        //site, moving dimension, direction, index
        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        int LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>::neigh_dir(int n, int d, bool dir, int mu) const {
            if (dir)
                return (n - n % (dim_mul_[d] * dimensions_[d]) +
                        (n + dim_mul_[d]) % (dim_mul_[d] * dimensions_[d])) * elem_per_site_ + mu;
            else
                return (n - n % (dim_mul_[d] * dimensions_[d]) +
                        (n - dim_mul_[d] + dim_mul_[d] * dimensions_[d]) % (dim_mul_[d] * dimensions_[d])) *
                       elem_per_site_ + mu;
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        void LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>::set_nearest_neighbours() {
            // int offset;
            for (uint i = 0; i < get_size(); i++) {
                // offset = 1;
                std::vector<T*> nn_of_site;
                // std::cout << "i: " << i << std::endl;
                for (uint d = 0; d < dimensions_.size(); d++) {
                    //std::cout << i-i%(offset*dimensions[d])+(i+offset)%(offset*dimensions[d]) << " - " << i-i%(offset*dimensions[d])+(i-offset+offset*dimensions[d])%(offset*dimensions[d]) << std::endl;
                    // x + nu
                    nn_of_site.push_back(&lattice_[neigh_dir(i, d, true,
                                                            0)]); // i-i%(offset*dimensions[d])+(i+offset)%(offset*dimensions[d])]);
                    // nn_of_site.push_back(&lattice_[i-i%(offset*dimensions[d])+(i+offset)%(offset*dimensions[d])]);
                    // x - nu
                    nn_of_site.push_back(&lattice_[neigh_dir(i, d, false, 0)]);
                    // nn_of_site.push_back(&lattice_[i-i%(offset*dimensions[d])+(i-offset+offset*dimensions[d])%(offset*dimensions[d])]);
                    // offset = offset*dimensions[d];
                }
                neighbours_.push_back(nn_of_site);
            }
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        void LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>::set_plaquette_neighbours() {
            // Loop over all sites
            for (uint n = 0; n < get_size() / elem_per_site_; n++) {
                // Loop over all links for a given site
                for (auto mu = 0; mu < dimension_; mu++) {
                    std::vector<T*> nn_of_site;
                    // Loop over possible plaquettes (left and right to the respective dimension)
                    for (auto nu = 0; nu < dimension_; nu++) {
                        if (nu != mu) {
                            /*lat.cumneigh_[lat.neigh_dir(n,mu,true,nu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                            lat.cumneigh_[n*lat.dim()+nu]++;
                            lat.cumneigh_[lat.neigh_dir(lat.neigh_dir(n,mu,true,0)/lat.dim(),nu,false,nu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,false,nu)]++;*/

                            nn_of_site.push_back(&lattice_[neigh_dir(n, mu, true, nu)]);
                            nn_of_site.push_back(&lattice_[neigh_dir(n, nu, true, mu)]);
                            nn_of_site.push_back(&lattice_[n * elem_per_site_ + nu]);
                            nn_of_site.push_back(
                                    &lattice_[neigh_dir(neigh_dir(n, mu, true, 0) / elem_per_site_, nu, false, nu)]);
                            nn_of_site.push_back(&lattice_[neigh_dir(n, nu, false, mu)]);
                            nn_of_site.push_back(&lattice_[neigh_dir(n, nu, false, nu)]);
                        }
                    }
                    neighbours_.push_back(nn_of_site);
                }
            }
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        std::ostream &operator<<(std::ostream &os,
                                 const LatticeSystem<T, Model, Method, UpdateDynamics, Sampler> &lattice) {
            for (uint i = 0; i < lattice.size(); i++) {
                std::cout << lattice(i) << " ";
                //if((i+1)%30 == 0) std::cout << std::endl;
            }
            std::cout << std::endl;
            return os;
        }

        template<typename T, typename Model, typename Method, typename UpdateDynamics, typename Sampler>
        std::vector<std::unique_ptr<mcmc::measures::Measure<LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>>>>
        LatticeSystem<T, Model, Method, UpdateDynamics, Sampler>::generate_lattice_system_measures(
                const std::vector<std::string> &measure_names) {
            typedef LatticeSystem<T, Model, Method, UpdateDynamics, Sampler> LatSys;
            std::vector<std::unique_ptr<mcmc::measures::Measure<LatSys>>> lattice_measures{};
            for (auto &measure_name :  measure_names)
                if (measure_name == "Energy")
                    lattice_measures.push_back(std::make_unique<util::system_measures::MeasureEnergyPolicy<LatSys>>());
                else if (measure_name == "Drift")
                    lattice_measures.push_back(std::make_unique<util::system_measures::MeasureDriftPolicy<LatSys>>());
            return lattice_measures;
        }

    }
}

#endif //MAIN_LATTICE_HPP
