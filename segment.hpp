#pragma once

#include "vector.hpp"
#include <algorithm>

namespace geom {

    template<size_t Dim, typename T>
    class Segment {
    public:
        using point_type = Point<Dim, T>;
		using vector_type = Vector<Dim, T>;
    private:
        point_type m_p1;
        point_type m_p2;

    public:
        Segment(const point_type& p1, const point_type& p2) : m_p1(p1), m_p2(p2) {
        }
        const point_type& p1() const { return m_p1; }
        const point_type& p2() const { return m_p2; }

        T length_sq() const {
            return (m_p2 - m_p1).length_sq();
        }

        T length() const {
            return (m_p2 - m_p1).length();
        }

        point_type project(const point_type& p) const {
            const vector_type v = m_p2 - m_p1;
            const T len_sq = v.length_sq();

            if (len_sq < Coord<T>::Epsilon) {
                return m_p1;
            }

            const T t = dot_product(p - m_p1, v) / len_sq;
            const T clamped_t = std::clamp(t, T(0), T(1));
            return m_p1 + v * clamped_t;
        }
    };

    template<size_t Dim, typename T>
    T distance(const Point<Dim, T>& p, const Segment<Dim, T>& s) {
        const Point<Dim, T> projected_point = s.project(p);
        return (p - projected_point).length();
    }

    template<size_t Dim, typename T>
    bool contains(const Point<Dim, T>& p, const Segment<Dim, T>& s) {
        return distance(p, s) < Coord<T>::Epsilon;
    }

    template<typename T> using Segment2 = Segment<2, T>;
    template<typename T> using Segment3 = Segment<3, T>;

    using Segment2d = Segment2<double>;
    using Segment3d = Segment3<double>;


} // namespace geom
