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
    // Вычисляем пропорции золотого сечения (alpha ≈ 0.618, beta ≈ 0.382)
    const double alpha = (std::sqrt(5.0) - 1.0) / 2.0; 
    const double beta = 1.0 - alpha; 
    
    // Начальное деление отрезка [a, b] на две пробные точки (Шаг 0)
    double l = a + beta * (b - a);   // Левая точка (lambda в конспекте)
    double m = a + alpha * (b - a);  // Правая точка (mu в конспекте)
    
    // Считаем значения функции в этих точках (вызываем тяжелую функцию f)
    double fl = f(l);
    double fm = f(m);

    // Основной цикл: пока длина отрезка больше заданной точности eps (Шаг 1)
    while ((b - a) > eps) {
        // Сравниваем значения функции, чтобы понять, какую часть отрезка отбросить
        if (fl > fm) {  
            // СЛУЧАЙ 1: Минимум точно лежит в правой части [l, b] (Шаг 2 в конспекте)
            a = l;      // Сдвигаем левую границу 'a' в точку 'l'
            l = m;      // Старая правая точка 'm' теперь становится левой точкой 'l'
            fl = fm;    // Значение в новой левой точке уже известно, просто копируем его
            
            // Считаем ТОЛЬКО одну новую правую точку
            m = a + alpha * (b - a); 
            fm = f(m);  // Единственный вызов функции на этой итерации
        } else {        
            // СЛУЧАЙ 2: Минимум точно лежит в левой части [a, m] (Шаг 3 в конспекте)
            b = m;      // Сдвигаем правую границу 'b' в точку 'm'
            m = l;      // Старая левая точка 'l' теперь становится правой точкой 'm'
            fm = fl;    // Значение в новой правой точке берем из старой левой
            
            // Считаем ТОЛЬКО одну новую левую точку
            l = a + beta * (b - a); 
            fl = f(l);  // Единственный вызов функции на этой итерации
        }
    }
    
    // Когда отрезок стал меньше eps, возвращаем его середину как ответ (Шаг 1, STOP)
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
        //====================== Step 1 ======================
        for (size_t j = 0; j < n; j++) {
            tau = golden_section_search([&](double t) { return f(y + t * D[j]); }, -0.1, 0.1, eps); // оптимизация по координатным осям
            y = y + tau * D[j]; // оптимальные перемещения по координатным осям
        }
        new_x = y; // y_(n+1)
        // euclide norm
        double residual = norm(new_x - x, 2); // residual

        if (residual <= eps) { // условие выхода
            x = new_x;
            break;
        }
        //====================== Step 2 ======================
        d = new_x - x; // x_(k+1) - x_k направление образца
        x = new_x; // x_(k+1)
        tau = golden_section_search([&](double t) { return f(x + t * d); }, 0, 3, eps); // оптимизация по направлению образца
        y = x + tau * d; // новое y_1
        iter++; // k = k + 1
    }
    return x;
}

} // namespace Opt
