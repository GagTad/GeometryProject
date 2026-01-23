#pragma once
#include "Triangulation.h"
#include <vector>
#include <set>
#include <map>
#include <cmath>

namespace Geom {

 
    template<typename T>
    struct NavTriangle {
        Triangle<T> tri;
        int id;
        std::vector<int> neighbors;  

        NavTriangle() : id(-1) {}
        NavTriangle(const Triangle<T>& t, int idx) : tri(t), id(idx) {}
    };

 
    template<typename T>
    struct Edge {
        Point<T> p1, p2;

        Edge(const Point<T>& a, const Point<T>& b) {
            
            if (a < b) {
                p1 = a;
                p2 = b;
            }
            else {
                p1 = b;
                p2 = a;
            }
        }

        bool operator==(const Edge& other) const {
            return (p1 == other.p1) && (p2 == other.p2);
        }

        bool operator<(const Edge& other) const {
            if (!(p1 == other.p1)) return p1 < other.p1;
            return p2 < other.p2;
        }
    };

    template<typename T>
    class NavMesh {
    private:
        std::vector<NavTriangle<T>> triangles;

        void buildAdjacency() {
 
            std::map<Edge<T>, std::vector<int>> edgeMap;

            for (size_t i = 0; i < triangles.size(); i++) {
                const Triangle<T>& tri = triangles[i].tri;

                Edge<T> e1(tri.a, tri.b);
                Edge<T> e2(tri.b, tri.c);
                Edge<T> e3(tri.c, tri.a);

                edgeMap[e1].push_back(i);
                edgeMap[e2].push_back(i);
                edgeMap[e3].push_back(i);
            }
 
            for (auto& pair : edgeMap) {
                const std::vector<int>& tris = pair.second;

                if (tris.size() == 2) {
                    int t1 = tris[0];
                    int t2 = tris[1];

                    triangles[t1].neighbors.push_back(t2);
                    triangles[t2].neighbors.push_back(t1);
                }
            }
        }

    public:
        NavMesh() {}

        void build(const std::vector<Triangle<T>>& tris) {
            triangles.clear();

            for (size_t i = 0; i < tris.size(); i++) {
                triangles.push_back(NavTriangle<T>(tris[i], i));
            }

            buildAdjacency();
        }

        int findTriangleContaining(const Point<T>& p) const {
            for (size_t i = 0; i < triangles.size(); i++) {
                if (triangles[i].tri.contains(p)) {
                    return i;
                }
            }
            return -1;
        }

        const NavTriangle<T>& getTriangle(int id) const {
            return triangles[id];
        }

        size_t size() const {
            return triangles.size();
        }

        const std::vector<NavTriangle<T>>& getTriangles() const {
            return triangles;
        }

        bool getSharedEdge(int t1Id, int t2Id, Point<T>& edgeStart, Point<T>& edgeEnd) const {
            if (t1Id < 0 || t1Id >= (int)triangles.size() ||
                t2Id < 0 || t2Id >= (int)triangles.size()) {
                return false;
            }

            const Triangle<T>& tri1 = triangles[t1Id].tri;
            const Triangle<T>& tri2 = triangles[t2Id].tri;

            Point<T> tri1Points[3] = { tri1.a, tri1.b, tri1.c };
            Point<T> tri2Points[3] = { tri2.a, tri2.b, tri2.c };

            for (int i = 0; i < 3; i++) {
                Point<T> p1 = tri1Points[i];
                Point<T> p2 = tri1Points[(i + 1) % 3];
                Edge<T> e1(p1, p2);

                for (int j = 0; j < 3; j++) {
                    Point<T> q1 = tri2Points[j];
                    Point<T> q2 = tri2Points[(j + 1) % 3];
                    Edge<T> e2(q1, q2);

                    if (e1 == e2) {
                        edgeStart = e1.p1;
                        edgeEnd = e1.p2;
                        return true;
                    }
                }
            }

            return false;
        }
    };

} // namespace Geom