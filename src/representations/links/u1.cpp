#include "../../../include/lattice_model_impl/representations/links/u1.hpp"


namespace lm_impl {
    namespace link {

        U1::U1(std::string init) {
            if (init == "null") {
                x_.push_back(std::complex<double>(0, 0));
            } else if (init == "identity") {
                x_.push_back(std::complex<double>(1, 0));
            } else {
                std::uniform_real_distribution<double> uniform(-1, 1);
                x_.push_back(std::complex<double>(uniform(mcmc::util::gen), uniform(mcmc::util::gen)));
                x_[0] = x_[0] / std::abs(x_[0]);
            }
        }

        U1::U1(double epsilon) {
            std::random_device rd;
            std::mt19937 gen{rd()};
            std::uniform_real_distribution<double> uniform(-epsilon, epsilon);

            x_.push_back(std::complex<double>(1, uniform(gen)));

            x_[0] = x_[0] / std::abs(x_[0]);

            /*double length;
            for(auto i = 1; i < 4; i++) {
                x_.push_back(uniform(gen));
                length += x_[i]*x_[i];
            }
            length = sqrt(length);
            for(int i = 1; i < 4; i++) x_[i] = epsilon*x_[i]/length;*/
        };


        U1::U1(std::complex<double> a) {
            x_.push_back(a);
        }


        std::complex<double> U1::trace() {
            return x_[0];
        }


        U1 &U1::adjungate() {
            x_[0] = std::conj(x_[0]);
            return *this;
        }

        U1 operator*(const U1 &x, const U1 &y) {
            return U1(x(0) * y(0));
        }

        U1 operator*(const U1 &x, const double &y) {
            return x(0) * y;
        }

        U1 operator-(const U1 &a, const U1 &b) {
            U1 temp(a);
            temp -= b;
            return temp;
        }

        std::ostream &operator<<(std::ostream &os, const U1 &x) {
            os << x(0);
            return os;
        }
    }
}

namespace std {
    std::string to_string(lm_impl::link::U1 x)
    {
        return std::to_string(x(0));
    }

    double fabs(lm_impl::link::U1 x)
    {
        return std::fabs(x(0));
    }
}