#ifndef MAIN_LATTICE_HPP
#define MAIN_LATTICE_HPP


#include <mcmc_simulation/header.hpp>
#include <param_helper/params.hpp>


#include "../util/measures/lattice_measures.hpp"


namespace lm_impl {
    namespace lattice_system {

        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        class LatticeSystem : public mcmc::simulation::MeasureSystemBase<LatticeSystem<T, Model, Method, UpdateDynamics>> {
        public:
            explicit LatticeSystem(const json params) : mcmc::simulation::MeasureSystemBase<LatticeSystem<T, Model, Method, UpdateDynamics>>(params)
            {
                dimensions = this->template get_entry<std::vector<int>>("dimensions", std::vector<int>{4, 4});
                lattice_action_type = this->template get_entry<std::string>("lattice_action_type", "nearest_neighbour");

                n_sites = 1;
                dim_mul.push_back(n_sites);
                for (auto dim: dimensions) {
                    n_sites *= dim;
                    dim_mul.push_back(n_sites);
                }
                dimension = dimensions.size();

                if(lattice_action_type == "nearest_neighbour")
                    elem_per_site = 1;
                else
                    elem_per_site = dimension;

                mcmc_model = Model(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics>, Model>(
                        *this, Model::name()));

                mcmc_method = Method(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics>, Method>(
                        *this, Method::name()));

                lattice_update = UpdateDynamics(
                    mcmc::util::generate_parameter_class_json<LatticeSystem<T, Model, Method, UpdateDynamics>, UpdateDynamics>(
                        *this, UpdateDynamics::name()));

                initialize_lattice();
                if (lattice_action_type == "plaquette")
                    set_plaquette_neighbours();
                else
                    set_nearest_neighbours();

                mcmc_model.init(*this);
                mcmc_method.init(*this);
                lattice_update.init(*this);
            }
            
            LatticeSystem(
                    Model mcmc_model_=Model(),
                    Method mcmc_method_=Method(),
                    UpdateDynamics lattice_update_=UpdateDynamics(),
                    const std::vector<int> dimensions_=std::vector<int> {4, 4},
                    const std::string lattice_action_type_="nearest_neighbour",
                    const std::vector<std::string> measures_=std::vector<std::string> {}
            ) : LatticeSystem(
                json {
                    {Model::name(), mcmc_model_.get_json()},
                    {Method::name(), mcmc_method_.get_json()},
                    {UpdateDynamics::name(), lattice_update_.get_json()},
                    {"dimensions", dimensions_},
                    {"lattice_action_type", lattice_action_type_},
                    {"measures", measures_}
                }
            )
            {}

            typedef LatticeSystem<T, Model, Method, UpdateDynamics> System;

            /* static std::string name() { -> definition not allowed since python program wont find the system parameters
                return "lattice";
            } */

            void write_to_file(const std::string rel_config_path) override {
                // ToDo: Write function in parameters that accepts various objects and does the same operation as is done here for an arbitrary number of parameters -> parameters can even be collected in a vector of std::vector<*Parameteres>

                std::string mcmc_model_params_path = this->template get_entry<std::string>(
                    Model::name() + "_path", rel_config_path);
                mcmc_model.write_to_file(mcmc_model_params_path);

                std::string mcmc_method_params_path = this->template get_entry<std::string>(
                    Method::name() + "_path", rel_config_path);
                mcmc_method.write_to_file(mcmc_method_params_path);

                std::string lattice_update_params_path = this->template get_entry<std::string>(
                    UpdateDynamics::name() + "_path", rel_config_path);
                lattice_update.write_to_file(lattice_update_params_path);

                json mcmc_model_parameters_ = mcmc_model.get_json();
                this->template delete_entry(Model::name());

                json update_parameters_ = mcmc_method.get_json();
                this->template delete_entry(Method::name());

                json lattice_update_parameters_ = lattice_update.get_json();
                this->template delete_entry(UpdateDynamics::name());

                param_helper::params::Parameters::write_to_file(rel_config_path, this->name());

                this->template add_entry(Model::name(), mcmc_model_parameters_);
                this->template add_entry(Method::name(), update_parameters_);
                this->template add_entry(UpdateDynamics::name(), lattice_update_parameters_);
            }

            param_helper::params::Parameters build_expanded_raw_parameters() const override {
                param_helper::params::Parameters parameters(this->params_);
                parameters.add_entry(Model::name(), mcmc_model.get_json());
                parameters.add_entry(Method::name(), mcmc_method.get_json());
                parameters.add_entry(UpdateDynamics::name(), lattice_update.get_json());
                return parameters;
            }

            const std::vector<int> &get_dimensions() const {
                return dimensions;
            }

            const auto &get_dimension() const {
                return dimension;
            }

            const uint &get_elems_per_site() const {
                return elem_per_site;
            }

            void initialize(std::string starting_mode) {
                // Needs to be called at the end so that generated objects (dynamics, mcmc_model, etc.) can already be used!
                this->generate_measures(this->measure_names());

                if(starting_mode == "hot")
                    for (auto &site : lattice)
                        site = mcmc_method.template random_state<T>();
                else // starting_mode == "cold"
                    for (auto &site : lattice)
                        site = mcmc_method.template cold_state<T>();
            }

            void update_step(uint measure_interval) {
                lattice_update(*this, measure_interval);
            }

            // Returns the total number of elements of the lattice - (equals the total number of sites if elem_per_site=1)
            auto get_size() const {
                return n_sites * elem_per_site;
            }

            auto at(int i) const {
                return lattice[i];
            }

            auto &at(int i) {
                return lattice[i];
            }

            auto get_system_representation() const {
                return lattice;
            }

            auto &get_system_representation() {
                return lattice;
            }

            auto get_elem_per_site() const {
                return elem_per_site;
            }

            auto &neighbours_at(int i) {
                return neighbours[i];
            }

            auto neighbours_at(int i) const {
                return neighbours[i];
            }

            auto &get_neighbours() {
                return neighbours;
            }

            auto get_neighbours() const {
                return neighbours;
            }

            void generate_measures(const std::vector<std::string> &measure_names) override {
                this->measurements_.clear();
                auto lattice_related_measures = generate_lattice_system_measures(this->measure_names());
                this->concat_measures(lattice_related_measures);

                auto mcmc_model_related_measures = mcmc_model. template generate_mcmc_model_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics>>(*this);
                this->concat_measures(mcmc_model_related_measures);

                auto mcmc_method_related_measures = mcmc_method. template generate_mcmc_method_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics>>(*this);
                this->concat_measures(mcmc_method_related_measures);

                auto lattice_update_related_measures = lattice_update. template generate_update_dynamics_measures<
                    LatticeSystem<T, Model, Method, UpdateDynamics>>(*this);
                this->concat_measures(lattice_update_related_measures);

                auto common_defined_measures = this->generate_systembase_measures(this->measure_names());
                this->concat_measures(common_defined_measures);
            }

            auto& get_mcmc_model() const {
                return mcmc_model;
            }

            auto& get_mcmc_model() {
                return mcmc_model;
            }

            auto& get_mcmc_method() const {
                return mcmc_method;
            }

            auto& get_mcmc_method() {
                return mcmc_method;
            }

            auto& get_lattice_update() const {
                return lattice_update;
            }

            auto& get_lattice_update() {
                return lattice_update;
            }

            auto energy() const {
                decltype(mcmc_model.get_energy_per_lattice_elem(lattice[0], neighbours[0])) energy(0);
                for (uint i = 0; i < get_size(); i++) {
                    energy += mcmc_model.get_energy_per_lattice_elem(lattice[i], neighbours[i]);
                }
                return energy;
            }

            auto drift_term() const {
                std::cerr << "Drift term computation needs to be implemented here" << std::endl;
                std::exit(EXIT_FAILURE);
                decltype(mcmc_model.get_potential(lattice[0], neighbours[0])) drift_term(0);
                for (uint i = 0; i < get_size(); i++) {
                    // ToDo: How can this be integrated?
                    // drift_term += mcmc_model.get_drift_term(lattice[i], neighbours[i]);
                }
                return drift_term;
            }

            void normalize(std::vector<T> &lattice_grid) {
                for (auto &elem : lattice_grid)
                    elem = mcmc_model.normalize(elem);
            }

            typedef T SiteType;

            // Only used for testing in simulation.hpp
            /* typedef Model MP_;
            typedef Method UP_; */

        protected:
            Model mcmc_model;
            Method mcmc_method;
            UpdateDynamics lattice_update;

            std::string lattice_action_type;

            uint n_sites; // Total number of sites
            uint elem_per_site; // Number of elements per site (is equal to dimension), but only for the link mcmc_model
            std::vector<int> dimensions; // Different dimensions
            std::vector<int> dim_mul; // Accumulated different dimensions (by product)
            int dimension; // Number of dimensions
            
            std::vector<T> lattice;
            std::vector<std::vector<T*> > neighbours;

            void initialize_lattice();

            int neigh_dir(int n, int d, bool dir, int mu) const;

            void set_nearest_neighbours();

            void set_plaquette_neighbours();

            std::vector<std::unique_ptr<mcmc::measures::Measure<LatticeSystem>>>
            generate_lattice_system_measures(const std::vector<std::string> &measure_names);
        };


        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        void LatticeSystem<T, Model, Method, UpdateDynamics>::initialize_lattice() {
            std::cout << "Lattice size: " << get_size() << std::endl;
            lattice = std::vector<T>(get_size(), T(0));
        }


        //site, moving dimension, direction, index
        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        int LatticeSystem<T, Model, Method, UpdateDynamics>::neigh_dir(int n, int d, bool dir, int mu) const {
            if (dir)
                return (n - n % (dim_mul[d] * dimensions[d]) +
                        (n + dim_mul[d]) % (dim_mul[d] * dimensions[d])) * elem_per_site + mu;
            else
                return (n - n % (dim_mul[d] * dimensions[d]) +
                        (n - dim_mul[d] + dim_mul[d] * dimensions[d]) % (dim_mul[d] * dimensions[d])) *
                       elem_per_site + mu;
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        void LatticeSystem<T, Model, Method, UpdateDynamics>::set_nearest_neighbours() {
            // int offset;
            for (uint i = 0; i < get_size(); i++) {
                // offset = 1;
                std::vector<T *> nn_of_site;
                // std::cout << "i: " << i << std::endl;
                for (uint d = 0; d < dimensions.size(); d++) {
                    //std::cout << i-i%(offset*dimensions[d])+(i+offset)%(offset*dimensions[d]) << " - " << i-i%(offset*dimensions[d])+(i-offset+offset*dimensions[d])%(offset*dimensions[d]) << std::endl;
                    // x + nu
                    nn_of_site.push_back(&lattice[neigh_dir(i, d, true,
                                                            0)]); // i-i%(offset*dimensions[d])+(i+offset)%(offset*dimensions[d])]);
                    // nn_of_site.push_back(&lattice[i-i%(offset*dimensions[d])+(i+offset)%(offset*dimensions[d])]);
                    // x - nu
                    nn_of_site.push_back(&lattice[neigh_dir(i, d, false, 0)]);
                    // nn_of_site.push_back(&lattice[i-i%(offset*dimensions[d])+(i-offset+offset*dimensions[d])%(offset*dimensions[d])]);
                    // offset = offset*dimensions[d];
                }
                neighbours.push_back(nn_of_site);
            }
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        void LatticeSystem<T, Model, Method, UpdateDynamics>::set_plaquette_neighbours() {
            // Loop over all sites
            for (uint n = 0; n < get_size() / elem_per_site; n++) {
                // Loop over all links for a given site
                for (auto mu = 0; mu < dimension; mu++) {
                    std::vector<T *> nn_of_site;
                    // Loop over possible plaquettes (left and right to the respective dimension)
                    for (auto nu = 0; nu < dimension; nu++) {
                        if (nu != mu) {
                            /*lat.cumneigh_[lat.neigh_dir(n,mu,true,nu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                            lat.cumneigh_[n*lat.dim()+nu]++;
                            lat.cumneigh_[lat.neigh_dir(lat.neigh_dir(n,mu,true,0)/lat.dim(),nu,false,nu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,false,nu)]++;*/

                            nn_of_site.push_back(&lattice[neigh_dir(n, mu, true, nu)]);
                            nn_of_site.push_back(&lattice[neigh_dir(n, nu, true, mu)]);
                            nn_of_site.push_back(&lattice[n * elem_per_site + nu]);
                            nn_of_site.push_back(
                                    &lattice[neigh_dir(neigh_dir(n, mu, true, 0) / elem_per_site, nu, false, nu)]);
                            nn_of_site.push_back(&lattice[neigh_dir(n, nu, false, mu)]);
                            nn_of_site.push_back(&lattice[neigh_dir(n, nu, false, nu)]);
                        }
                    }
                    neighbours.push_back(nn_of_site);
                }
            }
        }


        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        std::ostream &operator<<(std::ostream &os,
                                 const LatticeSystem<T, Model, Method, UpdateDynamics> &lattice) {
            for (uint i = 0; i < lattice.size(); i++) {
                std::cout << lattice(i) << " ";
                //if((i+1)%30 == 0) std::cout << std::endl;
            }
            std::cout << std::endl;
            return os;
        }

        template<typename T, typename Model, typename Method, typename UpdateDynamics>
        std::vector<std::unique_ptr<mcmc::measures::Measure<LatticeSystem<T, Model, Method, UpdateDynamics>>>>
        LatticeSystem<T, Model, Method, UpdateDynamics>::generate_lattice_system_measures(
                const std::vector<std::string> &measure_names) {
            typedef LatticeSystem<T, Model, Method, UpdateDynamics> LatSys;
            std::vector<std::unique_ptr<mcmc::measures::Measure<LatSys>>> lattice_measures{};
            for (auto &measure_name :  measure_names)
                if (measure_name == "Energy")
                    lattice_measures.push_back(std::make_unique<util::lattice_system_model_measures::MeasureEnergyPolicy<LatSys>>());
                else if (measure_name == "Drift")
                    lattice_measures.push_back(std::make_unique<util::lattice_system_model_measures::MeasureDriftPolicy<LatSys>>());
            return lattice_measures;
        }

    }
}

#endif //MAIN_LATTICE_HPP
