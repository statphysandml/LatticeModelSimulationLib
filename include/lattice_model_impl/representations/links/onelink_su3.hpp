#ifndef MAIN_OneLinkSU3_HPP
#define MAIN_OneLinkSU3_HPP

#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include <mcmc/mcmc_simulation/util/complex_type.hpp>

#include <lattice_model_impl/representations/link.hpp>
#include <lattice_model_impl/representations/links/su2.hpp>


using namespace std::literals;

namespace lm_impl {
    namespace link {

        double compute_x0(double beta, double a);

        class OneLinkSU3 : public Link<std::complex<double>> {
        private:
            using Link<std::complex<double>>::x_;
        public:
            OneLinkSU3(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d,
                        std::complex<double> e, std::complex<double> f, std::complex<double> g, std::complex<double> h,
                        std::complex<double> i);

            explicit OneLinkSU3(double epsilon);

            explicit OneLinkSU3(std::string init = "random");

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

        OneLinkSU3 operator*(const OneLinkSU3 &x, const OneLinkSU3 &y);

        OneLinkSU3 operator*(const OneLinkSU3 &x, const double &y);

        OneLinkSU3 operator-(const OneLinkSU3 &a, const OneLinkSU3 &b);

        OneLinkSU3 operator/(const OneLinkSU3 &x, const double &y);

        std::ostream &operator<<(std::ostream &os, const OneLinkSU3 &x);
    }
}


namespace std {
    std::string to_string(lm_impl::link::OneLinkSU3 x);

    double fabs(lm_impl::link::OneLinkSU3 x);

    std::array<std::complex<double>, 8> operator/(const std::array<std::complex<double>, 8>&x, const double &y);

    std::string to_string(const std::array<std::complex<double>, 8> x);
}

#endif //MAIN_SU3_HPP
