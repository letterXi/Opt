#ifndef MODIFIED_VECTOR_HPP
#define MODIFIED_VECTOR_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace MyMath {
template<typename T> class Vector : public std::vector<T> {
public:
    using std::vector<T>::vector;
    Vector operator+(const Vector& other) const {
        if (this->size() != other.size())
            throw std::logic_error("different sizes");
        Vector new_vector(this->size());
        for (size_t i = 0; i < this->size(); i++)
            new_vector[i] = other[i] + (*this)[i];
        return new_vector;
    }

    Vector operator-(const Vector& other) const {
        if (this->size() != other.size())
            throw std::logic_error("different sizes");
        Vector new_vector(this->size());
        for (size_t i = 0; i < this->size(); i++)
            new_vector[i] = (*this)[i] - other[i];
        return new_vector;
    }

    Vector operator*(T coeff) const {
        Vector new_vector(this->size());
        for (size_t i = 0; i < this->size(); i++)
            new_vector[i] = (*this)[i] * coeff;
        return new_vector;
    }
};

template<typename T> inline Vector<T> operator*(T coeff, const Vector<T>& v) {
    Vector<T> new_vector(v.size());
    for (size_t i = 0; i < v.size(); i++)
        new_vector[i] = v[i] * coeff;
    return new_vector;
}

template<typename T> inline std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); i++) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]";
    return os;
}

template<typename T> inline T norm(const Vector<T>& v, size_t n) {
    T res = 0;
    for (const auto& val : v)
        res += std::pow(val, n);
    return std::pow(res, 1.0 / static_cast<double>(n));
}
} // namespace MyMath
#endif
