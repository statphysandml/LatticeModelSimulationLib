//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP


#include "boost/numeric/odeint/stepper/symplectic_rkn_sb3a_mclachlan.hpp"
#include "boost/numeric/odeint/integrate/integrate_n_steps.hpp"

#include "../../mcmc_method_base.hpp"


namespace lm_impl {
    namespace mcmc_method {

        template<typename T, typename Model, typename Sampler>
        class HybridMonteCarloUpdate : public MCMCMethodBase<HybridMonteCarloUpdate<T, Model, Sampler>> {
        public:
            explicit HybridMonteCarloUpdate(const json params):
                MCMCMethodBase<HybridMonteCarloUpdate<T, Model, Sampler>>(params),
                dt_(this->template get_entry<double>("dt", 0.01)),
                n_(this->template get_entry<int>("n", 20)),
                m_(this->template get_entry<double>("m", 1.0))
            {
                normal_ = std::normal_distribution<double>(0.0, 1.0);
                rand_ = std::uniform_real_distribution<double>(0.0, 1.0);        
            }

            explicit HybridMonteCarloUpdate(const double dt=0.01, const int n=20, const double m=1.0
            ) : HybridMonteCarloUpdate(json{
                    {"dt", dt},
                    {"n",  n},
                    {"m", m}})
            {}

            template<typename System>
            void initialize(const System &system) {
                model_ptr_ = &system.get_mcmc_model();
                momenta_ = std::vector<double>(system.size(), 0.0);
                backup_momenta_ = std::vector<double>(system.size(), 0.0);
            }

            template<typename System>
            void operator()(System &system) {
                auto current_energy = system.energy();
                auto current_system_grid(system.get_system_representation());

                // Sample momenta
                std::generate(momenta_.begin(), momenta_.end(), [this]() { return normal_(mcmc::util::g_gen); });
                std::copy(momenta_.begin(), momenta_.end(), backup_momenta_.begin());

                HamiltonianSystemMomentum hamiltonian_system_momentum(model_ptr_, system.get_neighbours());
                HamiltonianSystemCoor hamiltonian_system_coor(m_);

                boost::numeric::odeint::integrate_n_steps(symplectic_stepper(),
                                                          std::make_pair(hamiltonian_system_coor, hamiltonian_system_momentum),
                                                          std::make_pair(boost::ref(system.get_system_representation()),
                                                                    boost::ref(momenta_)),
                                                          0.0, dt_, n_);

                system.normalize(system.get_system_representation());

                auto proposal_energy = system.energy();

                auto kinetic_term = std::inner_product(backup_momenta_.begin(), backup_momenta_.end(), backup_momenta_.begin(), 0.0);
                auto proposal_kinetic_term = std::inner_product(momenta_.begin(), momenta_.end(), momenta_.begin(), 0.0);

                // Accept/Reject step
                if (rand_(mcmc::util::g_gen) >= std::min(1.0, std::exp(
                        -1.0 * (proposal_energy - current_energy) - 0.5 * (proposal_kinetic_term - kinetic_term) / m_))) {
                    auto &system_grid = system.get_system_representation();
                    system_grid = current_system_grid;
                }
            }

            typedef boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan<std::vector<T> > symplectic_stepper;

        protected:

            struct HamiltonianSystemMomentum {
                HamiltonianSystemMomentum(const Model* model_ptr,
                    const std::vector<std::vector<T*>> &neighbours):
                    hmc_model_ptr_(model_ptr), hmc_neighbours_(neighbours) {}

                void operator()(const std::vector<T> &q, std::vector<T> &dpdt) const {
                    for (uint i = 0; i < q.size(); i++) {
                        dpdt[i] = -1.0 * hmc_model_ptr_->get_drift_term(q[i], hmc_neighbours_[i]);
                    }
                }

                const Model* hmc_model_ptr_;
                const std::vector<std::vector<T*>> &hmc_neighbours_;
            };

            struct HamiltonianSystemCoor {
                HamiltonianSystemCoor(const double m) : m_(m) {}

                void operator()(const std::vector<T> &p, std::vector<T> &dqdt) const {
                    for (uint i = 0; i < p.size(); i++) {
                        dqdt[i] = p[i] / m_;
                    }
                }

                const double m_;
            };


            double dt_;
            int n_;
            double m_;

            const Model* model_ptr_;

            std::vector<double> momenta_;
            std::vector<double> backup_momenta_;
            std::normal_distribution<double> normal_;
            std::uniform_real_distribution<double> rand_;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
