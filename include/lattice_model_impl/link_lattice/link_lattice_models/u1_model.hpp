#ifndef MAIN_U1_MODEL_HPP
#define MAIN_U1_MODEL_HPP


#include <lattice_model_impl/link_lattice/link_lattice_model.hpp>


namespace lm_impl {
    namespace link_lattice_system {
        class U1Model : public LinkLatticeModel<U1Model> {
        public:
            explicit U1Model(const json params):
                LinkLatticeModel(params),
                beta_(get_entry<double>("beta", 0.5))
            {}

            explicit U1Model(double beta=0.5) : U1Model(json{
                    {"beta", beta}
            })
            {}

            static const std::string type() {
                return "U1Model";
            }

            static uint N() {
                return 1;
            };

            template<typename T, typename T2=double_t>
            T2 get_potential(const T link, const std::vector<T*> neighbours) const {
                T A = calc_A(neighbours);
                return beta_ / U1Model::N() * (neighbours.size() / 3.0 - U1Model::N() *
                    std::real((link * A).trace()));
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T link, const std::vector<T*> neighbours) const {
                T A = calc_A(neighbours, false);
                return beta_ / U1Model::N() * (neighbours.size() / (2.0 * 3.0) - U1Model::N() *
                    std::real((link * A).trace()));
            }

        private:
            double beta_;
        };
    }
}

#endif //MAIN_U1_MODEL_HPP
