#include <lattice_model_impl/representations/links/onelink_su3.hpp>


namespace lm_impl {
    namespace link {

        OneLinkSU3::OneLinkSU3(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d,
                        std::complex<double> e, std::complex<double> f, std::complex<double> g, std::complex<double> h,
                        std::complex<double> i) {
            x_.push_back(a);
            x_.push_back(b);
            x_.push_back(c);
            x_.push_back(d);
            x_.push_back(e);
            x_.push_back(f);
            x_.push_back(g);
            x_.push_back(h);
            x_.push_back(i);
        }

        OneLinkSU3::OneLinkSU3(double epsilon) {
            // ( proposal state according to equation (4.30))
            std::uniform_real_distribution<double> distribution(-0.5, 0.5);

            x_.push_back(sgn(distribution(mcmc::util::random::g_gen)) * sqrt(1 - epsilon * epsilon));

            double length = 0;
            for (auto i = 1; i < 9; i++) {
                x_.push_back(distribution(mcmc::util::random::g_gen));
                length += std::abs(x_[i])*std::abs(x_[i]);
            }
            length = sqrt(length);
            for (int i = 1; i < 9; i++) x_[i] = epsilon * x_[i] / length;
        };

        OneLinkSU3::OneLinkSU3(std::string init) {
            if (init == "null") {
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
            } else if (init == "identity") {
                x_.push_back(1);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(1);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(1);
            } else {
                /*std::uniform_real_distribution<double> distribution(-1, 1);
                double length;
                while (true) {
                    x_.clear();
                    length = 0;
                    for (int i = 0; i < 9; i++) {
                        x_.push_back(distribution(mcmc::util::gen));
                        length += std::abs(x_[i]) * std::abs(x_[i]);
                    }
                    length = sqrt(length);
                    if (length <= 1)
                        break;
                }
                for (int i = 0; i < 9; i++) x_[i] = x_[i] / length;*/

                x_.push_back(1);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(1);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(0);
                x_.push_back(1);
            }
        }

        std::complex<double> OneLinkSU3::x0() {
            return x_[0];
        }

        std::complex<double> OneLinkSU3::x1() {
            return x_[1];
        }

        std::complex<double> OneLinkSU3::x2() {
            return x_[2];
        }

        std::complex<double> OneLinkSU3::x3() {
            return x_[3];
        }

        std::complex<double> OneLinkSU3::x4() {
            return x_[4];
        }

        std::complex<double> OneLinkSU3::x5() {
            return x_[5];
        }

        std::complex<double> OneLinkSU3::x6() {
            return x_[6];
        }

        std::complex<double> OneLinkSU3::x7() {
            return x_[7];
        }

        std::complex<double> OneLinkSU3::x8() {
            return x_[8];
        }

        OneLinkSU3 operator*(const OneLinkSU3 &x, const OneLinkSU3 &y) {
            return OneLinkSU3(x(0) * y(0) + x(1) * y(3), //has to be changed
                                x(0) * y(1) + x(1) * y(3),
                                x(1) * y(0) + x(3) * y(2),
                                x(1) * y(1) + x(3) * y(3), 0.0,0.0,0.0,0.0,0.0);
        }

        OneLinkSU3 operator*(const OneLinkSU3 &x, const double &y) {
            OneLinkSU3 x_(x);
            for (auto i = 0; i < 9; i++)
                x_(i) *= y;
            return x_;
        }

        OneLinkSU3 operator-(const OneLinkSU3 &a, const OneLinkSU3 &b) {
            OneLinkSU3 temp(a);
            temp -= b;
            return temp;
        }

        OneLinkSU3 operator/(const OneLinkSU3 &x, const double &y) {
            OneLinkSU3 temp(x);
            temp /= y;
            return temp;
        }

        std::ostream &operator<<(std::ostream &os, const OneLinkSU3 &x) {
            os << "((" << x(0) << "+I*" << x(3) << "," << x(2) << "+I*" << x(1) << "),(" << -x(2) << "+I*" << x(1)
               << "," << x(0) << "+I*" << -x(3) << "))";
            return os;
        }
    }
}

namespace std {
    std::string to_string(lm_impl::link::OneLinkSU3 x)
    {
        std::string conf = "";
        for(auto j = 0; j < 9; j++)
            conf += std::to_string(x(j)) + " ";
        conf = conf.substr(0, conf.size() -1);
        return conf;
    }

    double fabs(lm_impl::link::OneLinkSU3 x)
    {
        std::cerr << "Fabs not implemented for SU3 so far." << std::endl;
        std::exit(EXIT_FAILURE);
        return std::fabs(x(0));
    }

    std::array<std::complex<double>, 8> operator/(const std::array<std::complex<double>, 8>&x, const double &y) {
        std::array<std::complex<double>, 8> temp_(x);
        for (auto i = 0; i < 8; i++)
            temp_[i] /= y;
        return temp_;
    }

    std::string to_string(const std::array<std::complex<double>, 8> x)
    {
        std::string conf = "";
        for(auto j = 0; j < 8; j++)
            conf += std::to_string(x[j]) + " ";
        conf = conf.substr(0, conf.size() -1);
        return conf;
    }
}