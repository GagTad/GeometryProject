#pragma once
#include "Point.h"

namespace Geom {

    template<typename T>
    struct Segment {
        Point<T> p1, p2;

        Segment() {}
        Segment(const Point<T>& point1, const Point<T>& point2) : p1(point1), p2(point2) {}

        double length() const {
            return magnitude(p2 - p1);
        }
    };

    template<typename T>
    int orientation(const Point<T>& p, const Point<T>& q, const Point<T>& r) {
        T val = crossProduct(q - p, r - p);
        if constexpr (std::is_floating_point<T>::value) {
            if (std::abs(val) < EPS) return 0;
        }
        else {
            if (val == 0) return 0;
        }
        return (val > 0) ? 1 : -1;
    }

    template <typename T>
    bool on_segment(const Point<T>& p, const Point<T>& q, const Point<T>& r) {
        return (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
            q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y));
    }

    template <typename T>
    bool do_intersect(const Segment<T>& s1, const Segment<T>& s2) {
        Point<T> p1 = s1.p1, q1 = s1.p2;
        Point<T> p2 = s2.p1, q2 = s2.p2;

        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);

        if (o1 != o2 && o3 != o4) return true;

        if (o1 == 0 && on_segment(p1, p2, q1)) return true;
        if (o2 == 0 && on_segment(p1, q2, q1)) return true;
        if (o3 == 0 && on_segment(p2, p1, q2)) return true;
        if (o4 == 0 && on_segment(p2, q1, q2)) return true;

        return false;
    }

    template<typename T>
    Point<double> getIntersectionPoint(const Segment<T>& s1, const Segment<T>& s2) {
        double a1 = s1.p2.y - s1.p1.y;
        double b1 = s1.p1.x - s1.p2.x;
        double c1 = a1 * s1.p1.x + b1 * s1.p1.y;

        double a2 = s2.p2.y - s2.p1.y;
        double b2 = s2.p1.x - s2.p2.x;
        double c2 = a2 * s2.p1.x + b2 * s2.p1.y;

        double det = a1 * b2 - a2 * b1;
        if (std::abs(det) < EPS) return Point<double>(NAN, NAN); 

        return Point<double>((b2 * c1 - b1 * c2) / det, (a1 * c2 - a2 * c1) / det);
    }
}