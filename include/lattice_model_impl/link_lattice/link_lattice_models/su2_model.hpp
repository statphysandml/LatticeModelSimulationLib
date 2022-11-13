#ifndef MAIN_SU2_MODEL_HPP
#define MAIN_SU2_MODEL_HPP


#include <lattice_model_impl/link_lattice/link_lattice_model.hpp>


namespace lm_impl {
    namespace link_lattice_system {
        class SU2Model : public LinkLatticeModel<SU2Model> {
        public:
            explicit SU2Model(const json params):
                LinkLatticeModel(params),
                beta_(get_entry<double>("beta", 0.5))
            {}

            explicit SU2Model(double beta=0.5) : SU2Model(json{
                    {"beta", beta}
            })
            {}

            static const std::string type() {
                return "SU2Model";
            }

            static uint N() {
                return 2;
            };

            // Corresponds to the contribution of the link to the action: S[U] = \beta/N tr[neighbours.size() / 3 * identity - U * A]
            // See Gattringer - Quantum Chromodynamics on the Lattice Section 4.1.4
            template<typename T, typename T2=double_t>
            T2 get_potential(const T link, const std::vector<T*> neighbours) const {
                T A = calc_A(neighbours);
                return beta_ / SU2Model::N() * (neighbours.size() / 3.0 * SU2Model::N() -
                    std::real((link * A).trace()));
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T link, const std::vector<T*> neighbours) const {
                T A = calc_A(neighbours, false);
                return beta_ / SU2Model::N() * (neighbours.size() / (2.0 * 3.0) * SU2Model::N() -
                    std::real((link * A).trace()));
            }

        private:
            double beta_;
        };
    }
}

#endif //MAIN_SU2_MODEL_HPP
