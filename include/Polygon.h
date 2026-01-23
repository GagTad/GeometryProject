#pragma once
#include <vector>
#include <algorithm>

namespace Geom {

 
    template<typename T>
    struct AABB {
        T minX, minY, maxX, maxY;

        bool intersects(const AABB& other) const {
            return (minX <= other.maxX && maxX >= other.minX) &&
                (minY <= other.maxY && maxY >= other.minY);
        }
    };

    template<typename T>
    struct Polygon {
        std::vector<Point<T>> vertices;

        Polygon() {}
        Polygon(const std::vector<Point<T>>& verts) : vertices(verts) {}

        void addVertex(const Point<T>& p) {
            vertices.push_back(p);
        }
        size_t size() const { return vertices.size(); }

        AABB<T> getBounds() const {
            if (vertices.empty()) return { 0, 0, 0, 0 };

            T minX = vertices[0].x;
            T maxX = vertices[0].x;
            T minY = vertices[0].y;
            T maxY = vertices[0].y;

            for (const auto& v : vertices) {
                if (v.x < minX) minX = v.x;
                if (v.x > maxX) maxX = v.x;
                if (v.y < minY) minY = v.y;
                if (v.y > maxY) maxY = v.y;
            }
            return { minX, minY, maxX, maxY };
        }
    };

    template<typename T>
    bool isPointInPolygon(const Polygon<T>& poly, const Point<T>& p) {
        bool inside = false;
        size_t n = poly.vertices.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            const Point<T>& vi = poly.vertices[i];
            const Point<T>& vj = poly.vertices[j];

            if (((vi.y > p.y) != (vj.y > p.y)) &&
                (p.x < (vj.x - vi.x) * (p.y - vi.y) / (double)(vj.y - vi.y) + vi.x)) {
                inside = !inside;
            }
        }
        return inside;
    }
}