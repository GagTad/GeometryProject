#pragma once
#include "geometry.h"
#include <vector>
#include <algorithm>

namespace Geom {

    template<typename T>
    struct Triangle {
        Point<T> a, b, c;

        Triangle() {}
        Triangle(const Point<T>& p1, const Point<T>& p2, const Point<T>& p3)
            : a(p1), b(p2), c(p3) {
        }

        Point<T> centroid() const {
            return Point<T>((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
        }

        AABB<T> getBounds() const {
            return {
                std::min({a.x, b.x, c.x}),
                std::min({a.y, b.y, c.y}),
                std::max({a.x, b.x, c.x}),
                std::max({a.y, b.y, c.y})
            };
        }

        bool contains(const Point<T>& p) const {
            Vector<T> v0 = c - a;
            Vector<T> v1 = b - a;
            Vector<T> v2 = p - a;

            T dot00 = dotProduct(v0, v0);
            T dot01 = dotProduct(v0, v1);
            T dot02 = dotProduct(v0, v2);
            T dot11 = dotProduct(v1, v1);
            T dot12 = dotProduct(v1, v2);

            T invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
            T u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            T v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            return (u >= -EPS) && (v >= -EPS) && (u + v <= 1 + EPS);
        }
    };

    template<typename T>
    bool triangleIntersectsObstacle(const Triangle<T>& tri, const Polygon<T>& obstacle) {

        if (isPointInPolygon(obstacle, tri.centroid())) return true;
        if (isPointInPolygon(obstacle, tri.a)) return true;
        if (isPointInPolygon(obstacle, tri.b)) return true;
        if (isPointInPolygon(obstacle, tri.c)) return true;

        for (const auto& v : obstacle.vertices) {
            if (tri.contains(v)) return true;
        }
        //0.25 , 0,5 , 0,75 maseri bajanelov stugumner 
        Point<T> samples[9] = {
            Point<T>((tri.a.x + tri.b.x) / 2.0, (tri.a.y + tri.b.y) / 2.0),
            Point<T>((tri.a.x * 0.75 + tri.b.x * 0.25), (tri.a.y * 0.75 + tri.b.y * 0.25)),
            Point<T>((tri.a.x * 0.25 + tri.b.x * 0.75), (tri.a.y * 0.25 + tri.b.y * 0.75)),
            Point<T>((tri.b.x + tri.c.x) / 2.0, (tri.b.y + tri.c.y) / 2.0),
            Point<T>((tri.b.x * 0.75 + tri.c.x * 0.25), (tri.b.y * 0.75 + tri.c.y * 0.25)),
            Point<T>((tri.b.x * 0.25 + tri.c.x * 0.75), (tri.b.y * 0.25 + tri.c.y * 0.75)),
            Point<T>((tri.c.x + tri.a.x) / 2.0, (tri.c.y + tri.a.y) / 2.0),
            Point<T>((tri.c.x * 0.75 + tri.a.x * 0.25), (tri.c.y * 0.75 + tri.a.y * 0.25)),
            Point<T>((tri.c.x * 0.25 + tri.a.x * 0.75), (tri.c.y * 0.25 + tri.a.y * 0.75))
        };

        for (int i = 0; i < 9; i++) {
            if (isPointInPolygon(obstacle, samples[i])) return true;
        }
        //erankyan bolor koxmery argelqi bolor koxmeri het 
        Segment<T> triEdges[3] = {
            Segment<T>(tri.a, tri.b),
            Segment<T>(tri.b, tri.c),
            Segment<T>(tri.c, tri.a)
        };

        for (size_t i = 0; i < obstacle.size(); i++) {
            size_t next = (i + 1) % obstacle.size();
            Segment<T> obsEdge(obstacle.vertices[i], obstacle.vertices[next]);

            for (int j = 0; j < 3; j++) {
                if (do_intersect(triEdges[j], obsEdge)) {
                    bool sharedVertex =
                        (triEdges[j].p1 == obsEdge.p1 || triEdges[j].p1 == obsEdge.p2 ||
                            triEdges[j].p2 == obsEdge.p1 || triEdges[j].p2 == obsEdge.p2);

                    if (!sharedVertex) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    template<typename T>
    std::vector<Triangle<T>> triangulateWithHoles( const Polygon<T>& bounds, const std::vector<Polygon<T>>& obstacles) {

        std::vector<Triangle<T>> walkable;

        if (bounds.size() < 3) return walkable;

        T minX = bounds.vertices[0].x;
        T maxX = bounds.vertices[0].x;
        T minY = bounds.vertices[0].y;
        T maxY = bounds.vertices[0].y;

        for (const auto& v : bounds.vertices) {
            minX = std::min(minX, v.x);
            maxX = std::max(maxX, v.x);
            minY = std::min(minY, v.y);
            maxY = std::max(maxY, v.y);
        }

        std::vector<AABB<T>> obstacleBounds;
        for (const auto& obs : obstacles) {
            obstacleBounds.push_back(obs.getBounds());
        }


        int gridSize = 40;
        T dx = (maxX - minX) / gridSize;
        T dy = (maxY - minY) / gridSize;

        std::vector<std::vector<Point<T>>> grid;
        for (int j = 0; j <= gridSize; j++) {
            std::vector<Point<T>> row;
            for (int i = 0; i <= gridSize; i++) {
                T x = minX + i * dx;
                T y = minY + j * dy;
                row.push_back(Point<T>(x, y));
            }
            grid.push_back(row);
        }

        for (int j = 0; j < gridSize; j++) {
            for (int i = 0; i < gridSize; i++) {
                Point<T> p00 = grid[j][i];
                Point<T> p10 = grid[j][i + 1];
                Point<T> p01 = grid[j + 1][i];
                Point<T> p11 = grid[j + 1][i + 1];

                Triangle<T> tris[2] = {
                    Triangle<T>(p00, p10, p11),
                    Triangle<T>(p00, p11, p01)
                };

                for (int t = 0; t < 2; t++) {
                    const Triangle<T>& currentTri = tris[t];
                    Point<T> center = currentTri.centroid();

                    if (!isPointInPolygon(bounds, center)) continue;

                    bool valid = true;
                    AABB<T> triBox = currentTri.getBounds();

                    for (size_t k = 0; k < obstacles.size(); k++) {

                        if (!triBox.intersects(obstacleBounds[k])) {
                            continue;
                        }

                        if (triangleIntersectsObstacle(currentTri, obstacles[k])) {
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        walkable.push_back(currentTri);
                    }
                }
            }
        }

        return walkable;
    }

} // namespace Geom