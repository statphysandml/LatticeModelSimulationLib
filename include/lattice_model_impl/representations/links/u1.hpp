//
// Created by lukas on 04.11.19.
//

#ifndef MAIN_U1_HPP
#define MAIN_U1_HPP

#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include "mcmc_simulation/util/complex_type.hpp"
#include "../link.hpp"


namespace lm_impl {
    namespace link {

        class U1 : public Link<std::complex<double> > {
        private:
            using Link<std::complex<double> >::x_;
        public:
            U1(std::complex<double> a);

            U1(double epsilon);

            U1(std::string init = "random");

            std::complex<double> trace();

            U1 &adjungate() override;
        };

        U1 operator*(const U1 &x, const U1 &y);

        U1 operator*(const U1 &x, const double &y);

        template<typename T>
        U1 operator/(const U1 &x, const T &y) {
            U1 temp(x);
            temp /= y;
            return temp;
        }

        U1 operator-(const U1 &a, const U1 &b);

        std::ostream &operator<<(std::ostream &os, const U1 &x);
    }
}


namespace std {
    std::string to_string(lm_impl::link::U1 x);

    double fabs(lm_impl::link::U1 x);
}


#endif //MAIN_U1_HPP
