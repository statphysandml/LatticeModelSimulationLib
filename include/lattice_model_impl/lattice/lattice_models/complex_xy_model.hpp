#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP


#include <lattice_model_impl/lattice/mcmc_model_base.hpp>
#include <mcmc/mcmc_simulation/util/random.hpp>
#include <nlohmann/json.hpp>


namespace lm_impl {
    namespace lattice_system {

        class ComplexXYModel : public lm_impl::model::MCMCModelBase<ComplexXYModel> {
        public:
            explicit ComplexXYModel(const json params):
                MCMCModelBase(params),
                beta_(get_entry<double>("beta", 0.4)),
                mu_(get_entry<double>("mu", 0.0))
            {}

            explicit ComplexXYModel(double beta=0.4, double mu=0.0): ComplexXYModel(
                json{
                    {"beta", beta},
                    {"mu",   mu}
            }) {}

            static const std::string type() {
                return "ComplexXYModel";
            }

            std::complex<double> normalize(std::complex<double> state) {
                state.real(normalize(state.real()));
                return state;
            }

            double normalize(double state) {
                state = std::fmod(state, 2 * M_PI);
                if (state < 0) {
                    state += 2 * M_PI;
                }
                return state;
            }

            std::complex<double>
            get_potential(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) const
            {
                double S_re = 0;
                double S_im = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S_re += std::cos(site.real() - neighbours[i]->real()) *
                            std::cosh(site.imag() - neighbours[i]->imag() - mu_ * int(i == 0)) +
                            std::cos(neighbours[i + 1]->real() - site.real()) *
                            std::cosh(neighbours[i + 1]->imag() - site.imag() - mu_ * int(i == 0));
                    S_im += std::sin(site.real() - neighbours[i]->real()) *
                            std::sinh(site.imag() - neighbours[i]->imag() - mu_ * int(i == 0)) +
                            std::sin(neighbours[i + 1]->real() - site.real()) *
                            std::sinh(neighbours[i + 1]->imag() - site.imag() - mu_ * int(i == 0));
                }
                return {-1.0 * beta_ * S_re, beta_ * S_im};
            }

            std::complex<double>
            get_energy_per_lattice_elem(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) const
            {
                double S_re = 0;
                double S_im = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S_re += std::cos(site.real() - neighbours[i]->real()) *
                            std::cosh(site.imag() - neighbours[i]->imag() - mu_ * int(i == 0));
                    S_im += std::sin(site.real() - neighbours[i]->real()) *
                            std::sinh(site.imag() - neighbours[i]->imag() - mu_ * int(i == 0));
                }
                return {-1.0 * beta_ * S_re, beta_ * S_im};
            }

            std::complex<double>
            get_drift_term(const std::complex<double> site, const std::vector<std::complex<double> *> neighbours) const
            {
                double S_re = 0;
                double S_im = 0;
                for (size_t i = 0; i < neighbours.size(); i += 2) {
                    S_re += std::sin(site.real() - neighbours[i]->real()) *
                            std::cosh(site.imag() - neighbours[i]->imag() - mu_ * int(i == 0)) +
                            std::sin(site.real() - neighbours[i + 1]->real()) *
                            std::cosh(site.imag() - neighbours[i + 1]->imag() + mu_ * int(i == 0));
                    S_im += std::cos(site.real() - neighbours[i]->real()) *
                            std::sinh(site.imag() - neighbours[i]->imag() - mu_ * int(i == 0)) +
                            std::cos(site.real() - neighbours[i + 1]->real()) *
                            std::sinh(site.imag() - neighbours[i + 1]->imag() + mu_ * int(i == 0));
                }
                return beta_ * std::complex<double>(S_re, S_im);
            }

        private:
            double beta_;
            double mu_;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_XY_MODEL_HPP
