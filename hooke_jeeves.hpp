#pragma once
#include "my_math/vector.hpp"
#include <functional>
#include <iostream>

namespace Opt {
using Vector = MyMath::Vector<double>; // (x1, x2, x3, ..., xn)
using Vectors = MyMath::Vector<Vector>;
using FunctionND = std::function<double(const Vector&)>; // f(x1, ..., xn) -> R
using Function1D = std::function<double(double)>;        // f(x) -> R

// Golden section search
inline double golden_section_search(Function1D f, double a, double b, double eps) {
    const double alpha = (std::sqrt(5.0) - 1.0) / 2.0;
    const double beta = 1.0 - alpha;
    double l = a + beta * (b - a);
    double m = a + alpha * (b - a);

    double fl = f(l);
    double fm = f(m);

    while ((b - a) > eps) {
        if (fl > fm) {
            a = l;
            l = m;
            fl = fm;
            m = a + alpha * (b - a);
            fm = f(m);
        } else {
            b = m;
            m = l;
            fm = fl;
            l = a + beta * (b - a);
            fl = f(l);
        }
    }
    return (a + b) / 2.0;
}

// Create basis vectors
inline Vectors create_d(size_t n) {
    Vectors D(n);
    for (size_t i = 0; i < n; i++) {
        Vector temp(n, 0.0);
        temp[i] = 1.0;
        D[i] = temp;
    }
    return D;
}

inline Vector hooke_jeeves(FunctionND f, Vector x, double eps, bool verbose = false) {
    //============================================ Initial Stage ============================================
    size_t n = x.size();
    Vectors D = create_d(n); // e_j

    Vector new_x;
    Vector d;
    double tau;

    Vector y = x;
    size_t iter = 1;

    //============================================ Main Stage ============================================
    while (true) {
        x = y;
        if (verbose) {
            std::cout << "iteration " << iter << std::endl;
            std::cout << "x_" << iter << " = " << y << std::endl;
        }
        //====================== Step 1 ======================
        for (size_t j = 0; j < n; j++) {
            tau = golden_section_search([&](double t) { return f(y + t * D[j]); }, -0.1, 0.1, eps);
            y = y + tau * D[j];

            if (verbose)
                std::cout << " y_" << j+1 << " = " << y << std::endl;
        }
        new_x = y;
        // euclide norm
        double residual = norm(new_x - x, 2);

        if (verbose) {
            std::cout << "residual" << " = " << residual << std::endl << std::endl;
        }

        if (residual <= eps) {
            x = new_x;
            break;
        }
        //====================== Step 2 ======================
        d = new_x - x;
        x = new_x;
        tau = golden_section_search([&](double t) { return f(x + t * d); }, 0, 3, eps);
        y = x + tau * d;
        iter++;
    }
    if (verbose)
        std::cout << "x_min = " << x << std::endl;
    return x;
}

} // namespace Opt
