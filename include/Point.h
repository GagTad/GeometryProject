#pragma once
#include "Vector.h" 

namespace Geom {

    template<typename T>
    struct Point {
        T x, y;

        Point() : x(0), y(0) {}
        Point(T x_val, T y_val) : x(x_val), y(y_val) {}
    };

    template<typename T>
    Vector<T> operator-(const Point<T>& p1, const Point<T>& p2) {
        return Vector<T>(p1.x - p2.x, p1.y - p2.y);
    }

    template<typename T>
    Point<T> operator+(const Point<T>& p, const Vector<T>& v) {
        return Point<T>(p.x + v.dx, p.y + v.dy);
    }

    template<typename T>
    Point<T> operator-(const Point<T>& p, const Vector<T>& v) {
        return Point<T>(p.x - v.dx, p.y - v.dy);
    }


    template<typename T>
    bool operator==(const Point<T>& a, const Point<T>& b) {
        if constexpr (std::is_floating_point<T>::value) {
            return (std::abs(a.x - b.x) < EPS) && (std::abs(a.y - b.y) < EPS);
        }
        else {
            return a.x == b.x && a.y == b.y;
        }
    }

    template<typename T>
    bool operator<(const Point<T>& a, const Point<T>& b) {
        if constexpr (std::is_floating_point<T>::value) {
            if (std::abs(a.x - b.x) > EPS) return a.x < b.x;
            return a.y < b.y - EPS;
        }
        else {
            if (a.x != b.x) return a.x < b.x;
            return a.y < b.y;
        }
    }

    template<typename T>
    T distSq(const Point<T>& p1, const Point<T>& p2) {
        T dx = p1.x - p2.x;
        T dy = p1.y - p2.y;
        return dx * dx + dy * dy;
    }
}