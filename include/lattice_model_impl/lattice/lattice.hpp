//
// Created by lukas on 09.10.19.
//

#ifndef MAIN_LATTICE_HPP
#define MAIN_LATTICE_HPP

#include "mcmc_simulation/header.hpp"

#include "../util/measures/lattice_measures.hpp"

namespace lm_impl {
    namespace lattice_system {

        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        class LatticeSystem;

        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        class LatticeParameters : public mcmc::simulation::SystemBaseParameters {
        public:
            explicit LatticeParameters(const json params_) : SystemBaseParameters(params_) {

                dimensions = get_entry<std::vector<int> >("dimensions", std::vector<int>{4, 4});
                lattice_action_type = get_entry<std::string>("lattice_action_type", "nearest_neighbour");

                n_sites = 1;
                dim_mul.push_back(n_sites);
                for (auto dim: dimensions) {
                    n_sites *= dim;
                    dim_mul.push_back(n_sites);
                }
                dimension = dimensions.size();

                if (lattice_action_type == "nearest_neighbour")
                    elem_per_site = 1;
                else
                    elem_per_site = dimension;
                size = n_sites * elem_per_site;

                model_parameters = std::make_unique<ModelParameters>(
                        mcmc::util::generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, ModelParameters>(
                                *this, model_parameters->param_file_name()));

                update_parameters = std::make_unique<UpdateFormalismParameters>(
                        mcmc::util::generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, UpdateFormalismParameters>(
                                *this, update_parameters->param_file_name()));

                lattice_update_parameters = std::make_unique<LatticeUpdateFormalismParameters>(
                        mcmc::util::generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, LatticeUpdateFormalismParameters>(
                                *this, lattice_update_parameters->param_file_name()));
            }

            typedef LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> System;

            /* static std::string name() { -> definition not allowed since python program wont find the system parameters
                return "lattice";
            } */

            std::unique_ptr<System> generate() {
                return std::make_unique<System>(*this);
            }

            void write_to_file(const std::string &rel_config_path) override {
                // ToDo: Write function in parameters that accepts various objects and does the same operation as is done here for an arbitrary number of parameters -> parameters can even be collected in a vector of std::vector<*Parameteres>

                std::string model_params_path = get_entry<std::string>(model_parameters->param_file_name() + "_path",
                                                                       rel_config_path);
                model_parameters->write_to_file(model_params_path);

                std::string update_formalism_params_path = get_entry<std::string>(
                        update_parameters->param_file_name() + "_path", rel_config_path);
                update_parameters->write_to_file(update_formalism_params_path);

                std::string lattice_update_params_path = get_entry<std::string>(
                        lattice_update_parameters->param_file_name() + "_path", rel_config_path);
                lattice_update_parameters->write_to_file(lattice_update_params_path);

                json model_parameters_ = model_parameters->get_json();
                delete_entry(model_parameters->param_file_name());

                json update_parameters_ = update_parameters->get_json();
                delete_entry(update_parameters->param_file_name());

                json lattice_update_parameters_ = lattice_update_parameters->get_json();
                delete_entry(lattice_update_parameters->param_file_name());

                Parameters::write_to_file(rel_config_path, name());

                add_entry(model_parameters->param_file_name(), model_parameters_);
                add_entry(update_parameters->param_file_name(), update_parameters_);
                add_entry(lattice_update_parameters->param_file_name(), lattice_update_parameters_);
            }

            Parameters build_expanded_raw_parameters() const override {
                Parameters parameters(params);
                parameters.add_entry(model_parameters->param_file_name(), model_parameters->get_json());
                parameters.add_entry(update_parameters->param_file_name(), update_parameters->get_json());
                parameters.add_entry(lattice_update_parameters->param_file_name(),
                                     lattice_update_parameters->get_json());
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

            // Only used for testing in simulation.hpp
            typedef ModelParameters MP_;
            typedef UpdateFormalismParameters UP_;

        protected:
            template<typename, typename, typename, typename>
            friend
            class LatticeSystem;

            std::unique_ptr<ModelParameters> model_parameters;
            std::unique_ptr<UpdateFormalismParameters> update_parameters;
            std::unique_ptr<LatticeUpdateFormalismParameters> lattice_update_parameters;

            uint n_sites; // Total number of sites
            uint size; // Total number of elements on lattice
            std::vector<int> dimensions; // Different dimensions
            std::vector<int> dim_mul; // Accumulated different dimensions (by product)
            int dimension; // Number of dimensions
            uint elem_per_site; // Number of elements per site (is equal to dimension), but only for the link model

            std::string lattice_action_type;
        };


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        class LatticeSystem
                : public mcmc::simulation::SystemBase<LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> > {
        public:
            explicit LatticeSystem(
                    const LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lp_)
                    : lp(lp_) {
                model = std::make_unique<typename ModelParameters::Model>(*lp.model_parameters);
                update_formalism = std::make_unique<typename UpdateFormalismParameters::MCMCUpdate>(
                        *lp.update_parameters, *model);

                initialize_lattice();
                if (lp.lattice_action_type == "plaquette")
                    set_plaquette_neighbours();
                else
                    set_nearest_neighbours();

                lattice_update = std::make_unique<typename LatticeUpdateFormalismParameters::UpdateDynamics>(
                        *lp.lattice_update_parameters);

                update_formalism->initialize(*this);
                lattice_update->initialize(*this);
            }

            void update_step(uint measure_interval) {
                lattice_update->operator()(*this, measure_interval);
            }

            void initialize(std::string starting_mode) {
                // Needs to be called at the end so that generated objects (dynamics, model, etc.) can already be used!
                this->generate_measures(lp.measures);

                std::cout << "Note : Cold start not possible, so far" << std::endl;
            }

            // Returns the total number of elements of the lattice - (equals the total number of sites if elem_per_site=1)
            auto get_size() const {
                return lp.size;
            }

            auto &at(int i) {
                return lattice[i];
            }

            auto at(int i) const {
                return lattice[i];
            }

            auto get_system_representation() const {
                return lattice;
            }

            auto &get_system_representation() {
                return lattice;
            }

            auto get_elem_per_site() const {
                return lp.elem_per_site;
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

            void generate_measures(const json &measure_names) override {
                auto lattice_related_measures = generate_lattice_system_measures(lp.measures);
                this->concat_measures(lattice_related_measures);

                auto model_related_measures = model->template generate_model_measures<LatticeSystem,
                        LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>>(lp);
                this->concat_measures(model_related_measures);

                auto mcmc_update_related_measures = update_formalism->template generate_mcmc_update_measures<LatticeSystem,
                        LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>>(lp);
                this->concat_measures(mcmc_update_related_measures);

                auto lattice_update_related_measures = lattice_update->template generate_update_dynamics_measures<LatticeSystem,
                        LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>>(lp);
                this->concat_measures(lattice_update_related_measures);

                auto common_defined_measures = this->generate_systembase_measures(lp.measures);
                this->concat_measures(common_defined_measures);
            }

            auto &get_update_formalism() {
                return *update_formalism;
            }

            auto energy() const {
                decltype(model->get_energy_per_lattice_elem(lattice[0], neighbours[0])) energy(0);
                for (uint i = 0; i < get_size(); i++) {
                    energy += model->get_energy_per_lattice_elem(lattice[i], neighbours[i]);
                }
                return energy;
            }

            auto drift_term() const {
                std::cout << "Drift term computation needs to be implemented here" << std::endl;
                std::exit(EXIT_FAILURE);
                decltype(model->get_potential(lattice[0], neighbours[0])) drift_term(0);
                for (uint i = 0; i < get_size(); i++) {
                    // ToDo: How can this be integrated?
                    // drift_term += model->get_drift_term(lattice[i], neighbours[i]);
                }
                return drift_term;
            }

            void normalize(std::vector<T> &lattice_grid) {
                for (auto &elem : lattice_grid)
                    elem = model->normalize(elem);
            }

            typedef T SiteType;
        protected:
            const LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lp;

            std::vector<T> lattice;
            std::vector<std::vector<T *> > neighbours;

            std::unique_ptr<typename ModelParameters::Model> model;
            std::unique_ptr<typename UpdateFormalismParameters::MCMCUpdate> update_formalism;
            std::unique_ptr<typename LatticeUpdateFormalismParameters::UpdateDynamics> lattice_update;

            void initialize_lattice();

            int neigh_dir(int n, int d, bool dir, int mu) const;

            void set_nearest_neighbours();

            void set_plaquette_neighbours();

            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<LatticeSystem>>>
            generate_lattice_system_measures(const json &measure_names);
        };


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        void
        LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::initialize_lattice() {
            lattice = std::vector<T>(get_size(), T(0));
            for (auto &site : lattice)
                site = update_formalism->template random_state<T>();
        }


        //site, moving dimension, direction, index
        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        int
        LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::neigh_dir(int n,
                                                                                                                  int d,
                                                                                                                  bool dir,
                                                                                                                  int mu) const {
            if (dir)
                return (n - n % (lp.dim_mul[d] * lp.dimensions[d]) +
                        (n + lp.dim_mul[d]) % (lp.dim_mul[d] * lp.dimensions[d])) * lp.elem_per_site + mu;
            else
                return (n - n % (lp.dim_mul[d] * lp.dimensions[d]) +
                        (n - lp.dim_mul[d] + lp.dim_mul[d] * lp.dimensions[d]) % (lp.dim_mul[d] * lp.dimensions[d])) *
                       lp.elem_per_site + mu;
        }


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        void
        LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::set_nearest_neighbours() {
            // int offset;
            for (uint i = 0; i < get_size(); i++) {
                // offset = 1;
                std::vector<T *> nn_of_site;
                // std::cout << "i: " << i << std::endl;
                for (uint d = 0; d < lp.dimensions.size(); d++) {
                    //std::cout << i-i%(offset*lp.dimensions[d])+(i+offset)%(offset*lp.dimensions[d]) << " - " << i-i%(offset*lp.dimensions[d])+(i-offset+offset*lp.dimensions[d])%(offset*lp.dimensions[d]) << std::endl;
                    // x + nu
                    nn_of_site.push_back(&lattice[neigh_dir(i, d, true,
                                                            0)]); // i-i%(offset*lp.dimensions[d])+(i+offset)%(offset*lp.dimensions[d])]);
                    // nn_of_site.push_back(&lattice[i-i%(offset*lp.dimensions[d])+(i+offset)%(offset*lp.dimensions[d])]);
                    // x - nu
                    nn_of_site.push_back(&lattice[neigh_dir(i, d, false, 0)]);
                    // nn_of_site.push_back(&lattice[i-i%(offset*lp.dimensions[d])+(i-offset+offset*lp.dimensions[d])%(offset*lp.dimensions[d])]);
                    // offset = offset*lp.dimensions[d];
                }
                neighbours.push_back(nn_of_site);
            }
        }


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        void
        LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::set_plaquette_neighbours() {
            // Loop over all sites
            for (uint n = 0; n < get_size() / lp.elem_per_site; n++) {
                // Loop over all links for a given site
                for (auto mu = 0; mu < lp.dimension; mu++) {
                    std::vector<T *> nn_of_site;
                    // Loop over possible plaquettes (left and right to the respective dimension)
                    for (auto nu = 0; nu < lp.dimension; nu++) {
                        if (nu != mu) {
                            /*lat.cumneigh_[lat.neigh_dir(n,mu,true,nu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                            lat.cumneigh_[n*lat.dim()+nu]++;
                            lat.cumneigh_[lat.neigh_dir(lat.neigh_dir(n,mu,true,0)/lat.dim(),nu,false,nu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                            lat.cumneigh_[lat.neigh_dir(n,nu,false,nu)]++;*/

                            nn_of_site.push_back(&lattice[neigh_dir(n, mu, true, nu)]);
                            nn_of_site.push_back(&lattice[neigh_dir(n, nu, true, mu)]);
                            nn_of_site.push_back(&lattice[n * lp.elem_per_site + nu]);
                            nn_of_site.push_back(
                                    &lattice[neigh_dir(neigh_dir(n, mu, true, 0) / lp.elem_per_site, nu, false, nu)]);
                            nn_of_site.push_back(&lattice[neigh_dir(n, nu, false, mu)]);
                            nn_of_site.push_back(&lattice[neigh_dir(n, nu, false, nu)]);
                        }
                    }
                    neighbours.push_back(nn_of_site);
                }
            }
        }


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        std::ostream &operator<<(std::ostream &os,
                                 const LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lattice) {
            for (uint i = 0; i < lattice.size(); i++) {
                std::cout << lattice(i) << " ";
                //if((i+1)%30 == 0) std::cout << std::endl;
            }
            std::cout << std::endl;
            return os;
        }

        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
        std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>>>>
        LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::generate_lattice_system_measures(
                const json &measure_names) {
            typedef LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> LatSys;
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<LatSys>>> lattice_measures{};
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
