#pragma once
#include "Utils.h"

namespace Geom {

    template<typename T>
    struct Vector {
        T dx, dy;

        Vector() : dx(0), dy(0) {}
        Vector(T dx_val, T dy_val) : dx(dx_val), dy(dy_val) {}
    };

    // Vector operations
    template<typename T>
    Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2) {
        return Vector<T>(v1.dx + v2.dx, v1.dy + v2.dy);
    }

    template<typename T>
    Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2) {
        return Vector<T>(v1.dx - v2.dx, v1.dy - v2.dy);
    }

    template<typename T>
    Vector<T> operator*(const Vector<T>& v, T scalar) {
        return Vector<T>(v.dx * scalar, v.dy * scalar);
    }

    template<typename T>
    Vector<T> operator*(T scalar, const Vector<T>& v) {
        return Vector<T>(v.dx * scalar, v.dy * scalar);
    }

    template<typename T>
    Vector<T> operator/(const Vector<T>& v, T scalar) {
        return Vector<T>(v.dx / scalar, v.dy / scalar);
    }

    // Dot Product 
    template <typename T>
    T dotProduct(const Vector<T>& v1, const Vector<T>& v2) {
        return v1.dx * v2.dx + v1.dy * v2.dy;
    }

    // Cross Product 
    template <typename T>
    T crossProduct(const Vector<T>& v1, const Vector<T>& v2) {
        return v1.dx * v2.dy - v1.dy * v2.dx;
    }

    template<typename T>
    T magnitudeSq(const Vector<T>& v) {
        return v.dx * v.dx + v.dy * v.dy;
    }

    template<typename T >
    double magnitude(const Vector<T>& v) {
        return std::sqrt(static_cast<double>(magnitudeSq(v)));
    }


    template<typename T>
    double length(const Vector<T>& v) {
        return std::sqrt(v.dx * v.dx + v.dy * v.dy);
    }

    template<typename T>
    Vector<T> normalize(const Vector<T>& v) {
        double len = length(v);
        if (len < EPS) return Vector<T>(0, 0); 
        return Vector<T>(v.dx / len, v.dy / len);
    }
}