#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <numeric> 
#include <cassert>
#include "coord.hpp" 

namespace geom {

    template<size_t Dim, typename T>
    class Vector {
    public:
        using coord_type = Coord<T>;

    private:
        std::array<coord_type, Dim> m_coords; 

    public:
        Vector() {
            m_coords.fill(coord_type(0));
        }

   
        template<typename... Args>
        explicit Vector(Args... args) : m_coords{ coord_type(args)... } {
            static_assert(sizeof...(args) == Dim, "Incorrect number of arguments for Vector constructor.");
        }

        Vector(const std::initializer_list<T>& list) {
            assert(list.size() == Dim && "Incorrect number of arguments for initializer_list constructor.");
            size_t i = 0;
            for (const T& val : list) {
                m_coords[i++] = val;
            }
        }

        coord_type& operator[](size_t index) { return m_coords[index]; }
        const coord_type& operator[](size_t index) const { return m_coords[index]; }

        bool operator==(const Vector& other) const { return m_coords == other.m_coords; }
        bool operator!=(const Vector& other) const { return !(*this == other); }

        Vector& operator+=(const Vector& other) {
            for (size_t i = 0; i < Dim; ++i) { m_coords[i] += other.m_coords[i]; }
            return *this;
        }
        Vector& operator-=(const Vector& other) {
            for (size_t i = 0; i < Dim; ++i) { m_coords[i] -= other.m_coords[i]; }
            return *this;
        }

        Vector& operator*=(const T& scalar) {
            for (size_t i = 0; i < Dim; ++i) { m_coords[i] *= scalar; }
            return *this;
        }
        Vector& operator/=(const T& scalar) {
            for (size_t i = 0; i < Dim; ++i) { m_coords[i] /= scalar; }
            return *this;
        }

        T length_sq() const {
            T result = 0;
            for (size_t i = 0; i < Dim; ++i) {
                result += m_coords[i].value * m_coords[i].value;
            }
            return result;
        }

        T length() const {
            return std::sqrt(length_sq());
        }

        void normalize() {
            T len = length();
            if (len > 0) { 
                *this /= len;
            }
        }
    };

    template<size_t Dim, typename T>
    Vector<Dim, T> operator+(Vector<Dim, T> lhs, const Vector<Dim, T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<size_t Dim, typename T>
    Vector<Dim, T> operator-(Vector<Dim, T> lhs, const Vector<Dim, T>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<size_t Dim, typename T>
    Vector<Dim, T> operator*(Vector<Dim, T> vec, const T& scalar) {
        vec *= scalar;
        return vec;
    }
    template<size_t Dim, typename T>
    Vector<Dim, T> operator*(const T& scalar, Vector<Dim, T> vec) {
        vec *= scalar;
        return vec;
    }

    template<size_t Dim, typename T>
    T dot_product(const Vector<Dim, T>& a, const Vector<Dim, T>& b) {
        T result = 0;
        for (size_t i = 0; i < Dim; ++i) {
            result += a[i].value * b[i].value;
        }
        return result;
    }

    template<size_t Dim, typename T>
    std::ostream& operator<<(std::ostream& os, const Vector<Dim, T>& vec) {
        os << "Vector<" << Dim << ">(";
        for (size_t i = 0; i < Dim; ++i) {
            os << vec[i] << (i == Dim - 1 ? "" : ", ");
        }
        os << ")";
        return os;
    }

    template<size_t Dim, typename T> using Point = Vector<Dim, T>;

    template<typename T> using Vector2 = Vector<2, T>;
    template<typename T> using Vector3 = Vector<3, T>;

    template<typename T> using Point2 = Vector<2, T>;
    template<typename T> using Point3 = Vector<3, T>;

    using Vector2d = Vector2<double>;
    using Vector3d = Vector3<double>;
    using Point2d = Point2<double>;
    using Point3d = Point3<double>;

} // namespace geom
