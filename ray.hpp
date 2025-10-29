#pragma once

#include <cassert>
#include "vector.hpp"
#include "line.hpp"

namespace geom {

    template<size_t Dim, typename T>
    class Ray {
    public:
        using point_type = Point<Dim, T>;
        using vector_type = Vector<Dim, T>;

    private:
        point_type m_origin;
        vector_type m_direction;

        Ray(const point_type& origin, const vector_type& direction)
            : m_origin(origin), m_direction(direction) {
            this->m_direction.normalize();
        }

    public:
        static Ray from_point_direction(const point_type& origin, const vector_type& direction) {
            assert(direction.length_sq() > 0 && "Ray direction vector cannot be zero.");
            return Ray(origin, direction);
        }

        static Ray from_points(const point_type& start_point, const point_type& second_point) {
            assert((second_point - start_point).length_sq() > 0 && "Points for Ray construction cannot be the same.");
            return Ray(start_point, second_point - start_point);
        }

        const point_type& origin() const { return this->m_origin; }
        const vector_type& direction() const { return this->m_direction; }

        point_type point_at(const T& t) const {
            assert(t >= 0 && "Parameter t for Ray::point_at must be non-negative.");
            return this->m_origin + (this->m_direction * t);
        }

        bool contains(const point_type& p) const {
            Line<Dim, T> supporting_line = Line<Dim, T>::from_point_direction(m_origin, m_direction);
            if (!supporting_line.contains(p)) {
                return false;
            }
            return dot_product(p - m_origin, m_direction) >= 0;
        }
    };

    template<typename T> using Ray2 = Ray<2, T>;
    template<typename T> using Ray3 = Ray<3, T>;

    using Ray2d = Ray2<double>;
    using Ray3d = Ray3<double>;


} // namespace geom
