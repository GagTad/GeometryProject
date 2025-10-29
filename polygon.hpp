#pragma once

#include <vector>
#include <cassert>
#include "vector.hpp"
#include "segment.hpp"
#include "algorithms.hpp" // cross_product-ի համար

namespace geom {

    template<size_t Dim, typename T>
    class Polygon {
    public:
        using point_type = Point<Dim, T>;

    private:
        std::vector<point_type> m_vertices;

    public:
        // --- Կոնստրուկտոր ---
        Polygon(const std::vector<point_type>& vertices) : m_vertices(vertices) {
            assert(m_vertices.size() >= 3 && "Polygon must have at least 3 vertices.");
        }

        // --- Հասանելիության մեթոդներ ---
        size_t num_vertices() const {
            return m_vertices.size();
        }

        const std::vector<point_type>& vertices() const {
            return m_vertices;
        }

        // Վերադարձնում է i-րդ կողմը (հատվածը)
        Segment<Dim, T> edge(size_t i) const {
            assert(i < num_vertices() && "Edge index out of bounds.");
            // Վերջին կողմը միացնում է վերջին գագաթը առաջինին
            return Segment<Dim, T>(m_vertices[i], m_vertices[(i + 1) % num_vertices()]);
        }

        // --- Երկրաչափական մեթոդներ (Ալգորիթմներ) ---

        // 2D բազմանկյան մակերեսի հաշվարկ (Shoelace formula)
        T area() const {
            // Այս մեթոդը իմաստ ունի միայն 2D-ի համար
            static_assert(Dim == 2, "Area calculation is only implemented for 2D polygons.");

            T total_area = 0;
            for (size_t i = 0; i < num_vertices(); ++i) {
                point_type p1 = m_vertices[i];
                point_type p2 = m_vertices[(i + 1) % num_vertices()];
                total_area += cross_product(p1, p2);
            }
            return std::abs(total_area) / 2.0;
        }
    };

    // --- Տիպերի այլանուններ ---
    template<typename T> using Polygon2 = Polygon<2, T>;
    template<typename T> using Polygon3 = Polygon<3, T>;

    using Polygon2d = Polygon2<double>;
    using Polygon3d = Polygon3<double>;

} // namespace geom