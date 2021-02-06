//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP


#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../../mcmc_update_base.hpp"


namespace lm_impl {
    namespace mcmc_update {


        template<typename T, typename ModelParameters, typename SamplerCl>
        class HybridMonteCarloUpdate;


        template<typename T, typename ModelParameters, typename SamplerCl>
        class HybridMonteCarloUpdateParameters : public MCMCUpdateBaseParameters {
        public:
            explicit HybridMonteCarloUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                            dt(get_entry<double>("dt", 0.01)),
                                                                            n(get_entry<int>("n", 20)),
                                                                            m(get_entry<double>("m", 1.0)) {}

            explicit HybridMonteCarloUpdateParameters(const double dt_, const int n_, const double m_=1.0
            ) : HybridMonteCarloUpdateParameters(json{
                    {"dt", dt_},
                    {"n",  n_},
                    {"m", m_}}) {}

            static std::string name() {
                return "HybridMonteCarloUpdate";
            }

            typedef HybridMonteCarloUpdate<T, ModelParameters, SamplerCl> MCMCUpdate;
            typedef typename ModelParameters::Model Model;

        protected:
            friend class HybridMonteCarloUpdate<T, ModelParameters, SamplerCl>;

            const double dt;
            const int n;
            const double m;
        };


        template<typename T, typename ModelParameters, typename SamplerCl>
        class HybridMonteCarloUpdate
                : public MCMCUpdateBase<HybridMonteCarloUpdate<T, ModelParameters, SamplerCl>, SamplerCl> {
        public:
            explicit HybridMonteCarloUpdate(const HybridMonteCarloUpdateParameters<T, ModelParameters, SamplerCl> &up_,
                                            typename ModelParameters::Model &model_)
                    : MCMCUpdateBase<HybridMonteCarloUpdate<T, ModelParameters, SamplerCl>, SamplerCl>(up_.eps),
                      up(up_), model(model_) {
                normal = std::normal_distribution<double>(0.0, 1.0);
                rand = std::uniform_real_distribution<double>(0, 1);
            }

            template<typename Lattice>
            void initialize_mcmc_update(const Lattice &lattice) {
                momenta = std::vector<double>(lattice.size(), 0.0);
                backup_momenta = std::vector<double>(lattice.size(), 0.0);
            }

            template<typename Lattice>
            void operator()(Lattice &lattice) {
                auto current_energy = lattice.energy();
                auto current_lattice_grid(lattice.get_system_representation());

                // Sample momenta
                std::generate(momenta.begin(), momenta.end(), [this]() { return normal(mcmc::util::gen); });
                std::copy(momenta.begin(), momenta.end(), backup_momenta.begin());

                HamiltonianSystemMomentum hamiltonian_system_momentum(model, lattice.get_neighbours());
                HamiltonianSystemCoor hamiltonian_system_coor(up.m);

                boost::numeric::odeint::integrate_n_steps(symplectic_stepper(),
                                                          std::make_pair(hamiltonian_system_coor, hamiltonian_system_momentum),
                                                          std::make_pair(boost::ref(lattice.get_system_representation()),
                                                                    boost::ref(momenta)),
                                                          0.0, up.dt, up.n);

                lattice.normalize(lattice.get_system_representation());

                auto proposal_energy = lattice.energy();

                auto kinetic_term = std::inner_product(backup_momenta.begin(), backup_momenta.end(), backup_momenta.begin(), 0.0);
                auto proposal_kinetic_term = std::inner_product(momenta.begin(), momenta.end(), momenta.begin(), 0.0);

                // Accept/Reject step
                if (rand(mcmc::util::gen) >= std::min(1.0, std::exp(
                        -1.0 * (proposal_energy - current_energy) - 0.5 * (proposal_kinetic_term - kinetic_term) / up.m))) {
                    auto &lattice_grid = lattice.get_system_representation();
                    lattice_grid = current_lattice_grid;
                }
            }

            typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan<std::vector<T> > symplectic_stepper;

        protected:

            struct HamiltonianSystemMomentum {
                HamiltonianSystemMomentum(typename ModelParameters::Model &model_,
                                  const std::vector<std::vector<T *> > &neighbours_) : hmc_model(model_),
                                                                                       hmc_neighbours(neighbours_) {}

                void operator()(const std::vector<T> &q, std::vector<T> &dpdt) const {
                    for (uint i = 0; i < q.size(); i++) {
                        dpdt[i] = -1.0 * hmc_model.get_drift_term(q[i], hmc_neighbours[i]);
                    }
                }

                typename ModelParameters::Model &hmc_model;
                const std::vector<std::vector<T *> > &hmc_neighbours;
            };

            struct HamiltonianSystemCoor {
                HamiltonianSystemCoor(const double m_) : m(m_) {}

                void operator()(const std::vector<T> &p, std::vector<T> &dqdt) const {
                    for (uint i = 0; i < p.size(); i++) {
                        dqdt[i] = p[i] / m;
                    }
                }

                const double m;
            };

            const HybridMonteCarloUpdateParameters<T, ModelParameters, SamplerCl> &up;
            typename ModelParameters::Model &model;

            std::vector<double> momenta;
            std::vector<double> backup_momenta;
            std::normal_distribution<double> normal;
            std::uniform_real_distribution<double> rand;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
