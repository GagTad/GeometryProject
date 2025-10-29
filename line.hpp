#pragma once

#include <cassert>
#include "vector.hpp"

namespace geom {

    template<size_t Dim, typename T>
    class Line {
    public:
        using point_type = Point<Dim, T>;
        using vector_type = Vector<Dim, T>;

    private:
 
        Line(const point_type& origin, const vector_type& direction)
            : m_origin(origin), m_direction(direction) {
            this->m_direction.normalize();
        }

        point_type m_origin;
        vector_type m_direction;

    public:
        static Line from_points(const point_type& p1, const point_type& p2) {
            vector_type dir = p2 - p1;
            assert(dir.length_sq() > 0 && "Points for Line construction cannot be the same.");
            return Line(p1, dir); 
        }

        static Line from_point_direction(const point_type& origin, const vector_type& direction) {
            assert(direction.length_sq() > 0 && "Line direction vector cannot be zero.");
            return Line(origin, direction); 
        }

        const point_type& origin() const { return this->m_origin; }
        const vector_type& direction() const { return this->m_direction; }

        point_type point_at(const T& t) const {
            return this->m_origin + (this->m_direction * t);
        }

        point_type project(const point_type& p) const {
            vector_type to_p = p - this->m_origin;
            T t = dot_product(to_p, this->m_direction);
            return this->point_at(t);
        }

        bool contains(const point_type& p) const {
            return p == this->project(p);
        }
    };

    template<typename T> using Line2 = Line<2, T>;
    template<typename T> using Line3 = Line<3, T>;

    using Line2d = Line2<double>;
    using Line3d = Line3<double>;


} // namespace geom
