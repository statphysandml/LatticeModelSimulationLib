#ifndef MAIN_OneLinkSU3_HPP
#define MAIN_OneLinkSU3_HPP

#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include "../link.hpp"

using namespace std::literals;

namespace lm_impl {

    namespace link {

        double compute_x0(double beta, double a);

        template<typename T>
        class OneLinkSU3 : public Link<T> {
        private:
            using Link<T>::x_;
        public:
            OneLinkSU3(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d,
                        std::complex<double> e, std::complex<double> f, std::complex<double> g, std::complex<double> h,
                        std::complex<double> i);

            explicit OneLinkSU3(double epsilon);

            explicit OneLinkSU3(std::string init = "random");

            OneLinkSU3(const OneLinkSU3 &A, double beta);

            std::complex<double> trace();

            std::complex<double> det();

            OneLinkSU3 &adjungate();

            std::complex<double> x0();
            std::complex<double> x1();
            std::complex<double> x2();
            std::complex<double> x3();
            std::complex<double> x4();
            std::complex<double> x5();
            std::complex<double> x6();
            std::complex<double> x7();
            std::complex<double> x8();
        };

        template<typename T>
        OneLinkSU3<T>::OneLinkSU3(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d,
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


        template<typename T>
        OneLinkSU3<T>::OneLinkSU3(double epsilon) {
            // ( proposal state according to equation (4.30))
            std::uniform_real_distribution<double> distribution(-0.5, 0.5);

            x_.push_back(sgn(distribution(mcmc::util::g_gen)) * sqrt(1 - epsilon * epsilon));

            double length = 0;
            for (auto i = 1; i < 9; i++) {
                x_.push_back(distribution(mcmc::util::g_gen));
                length += std::abs(x_[i])*std::abs(x_[i]);
            }
            length = sqrt(length);
            for (int i = 1; i < 9; i++) x_[i] = epsilon * x_[i] / length;
        };

        template<typename T>
        OneLinkSU3<T>::OneLinkSU3(std::string init) {
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

        template<typename T>
        OneLinkSU3<T>::OneLinkSU3(const OneLinkSU3<T> &A, double beta) {
            // Proposal according to heatbath!!
            *this = A;
            double determinante = this->det();
            if (determinante == 0)
                *this = OneLinkSU3<T>();
            else {
                double a = sqrt(determinante);
                *this /= a; // = V
                double x0 = compute_x0(beta, a);
                double length_x = sqrt(1 - x0 * x0);
                std::uniform_real_distribution<double> distribution(-1, 1);
                double x1, x2, x3;
                double length_r;
                while (true) {
                    x1 = distribution(mcmc::util::g_gen);
                    x2 = distribution(mcmc::util::g_gen);
                    x3 = distribution(mcmc::util::g_gen);
                    length_r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                    if (length_r <= 1)
                        break;
                }
                SU2<T> X(x0, x1 * length_x / length_r, x2 * length_x / length_r, x3 * length_x / length_r);
                *this = X * (this->adjungate());
            }
        }

        template<typename T>                                // These functions are not adjusted to SU3!!
        std::complex<double> OneLinkSU3<T>::trace() {
            return x_[0] * 2.0;
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::det() {
            double determinante = 0;
            for (auto i = 0; i < 4; i++) determinante += x_[i] * x_[i];
            return determinante;
        }

        template<typename T>
        OneLinkSU3<T> &OneLinkSU3<T>::adjungate() {
            x_[3] = -x_[3];
            x_[2] = -x_[2];
            x_[1] = -x_[1];
            return *this;
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x0() {
            return x_[0];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x1() {
            return x_[1];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x2() {
            return x_[2];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x3() {
            return x_[3];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x4() {
            return x_[4];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x5() {
            return x_[5];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x6() {
            return x_[6];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x7() {
            return x_[7];
        }

        template<typename T>
        std::complex<double> OneLinkSU3<T>::x8() {
            return x_[8];
        }
        

        template<typename T>
        OneLinkSU3<T> operator*(const OneLinkSU3<T> &x, const double &y) {
            OneLinkSU3<T> x_(x);
            for (auto i = 0; i < 9; i++)
                x_(i) *= y;
            return x_;
        }

        template<typename T>
        OneLinkSU3<T> operator-(const OneLinkSU3<T> &a, const OneLinkSU3<T> &b) {
            OneLinkSU3<T> temp(a);
            temp -= b;
            return temp;
        }

        template<typename T, typename T2>
        OneLinkSU3<T> operator/(const OneLinkSU3<T> &x, const T2 &y) {
            OneLinkSU3<T> temp(x);
            temp /= y;
            return temp;
        }

        template<typename T>
        OneLinkSU3<T> operator*(const OneLinkSU3<T> &x, const OneLinkSU3<T> &y) {

            

            return OneLinkSU3<T>(x(0) * y(0) + x(1) * y(3), //has to be changed
                                x(0) * y(1) + x(1) * y(3),
                                x(1) * y(0) + x(3) * y(2),
                                x(1) * y(1) + x(3) * y(3), 0.0,0.0,0.0,0.0,0.0);
        }

        template<typename T>
        std::ostream &operator<<(std::ostream &os, const OneLinkSU3<T> &x) {
            os << "((" << x(0) << "+I*" << x(3) << "," << x(2) << "+I*" << x(1) << "),(" << -x(2) << "+I*" << x(1)
               << "," << x(0) << "+I*" << -x(3) << "))";
            return os;
        }

    }
}


namespace std {
    template<typename T>
    std::string to_string(lm_impl::link::OneLinkSU3<T> x)
    {
        std::string conf = "";
        for(auto j = 0; j < 9; j++)
            conf += std::to_string(x(j)) + " ";
        conf = conf.substr(0, conf.size() -1);
        return conf;
    }

    template<typename T>
    double fabs(lm_impl::link::OneLinkSU3<T> x)
    {
        std::cerr << "Fabs not implemented for SU3 so far." << std::endl;
        std::exit(EXIT_FAILURE);
        return std::fabs(x(0));
    }
}

#endif //MAIN_SU3_HPP
