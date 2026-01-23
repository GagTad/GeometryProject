#pragma once
#include "NavMesh.h"
#include <queue>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>

namespace Geom {

    struct AStarNode {
        int triangleId;
        double gCost; 
        double hCost;  
        int parent;

        AStarNode() : triangleId(-1), gCost(0), hCost(0), parent(-1) {}
        AStarNode(int id, double g, double h, int p)
            : triangleId(id), gCost(g), hCost(h), parent(p) {
        }

        double fCost() const { return gCost + hCost; }

        bool operator>(const AStarNode& other) const {
            return fCost() > other.fCost();
        }
    };

    template<typename T>
    class AStarPathfinder {
    private:
        const NavMesh<T>* navMesh;

        double heuristic(const Point<T>& p1, const Point<T>& p2) const {
            return std::sqrt(distSq(p1, p2));
        }

    public:
        AStarPathfinder(const NavMesh<T>* mesh) : navMesh(mesh) {}

        std::vector<int> findPath(const Point<T>& start, const Point<T>& goal) {
            std::vector<int> path;

            if (!navMesh || navMesh->size() == 0) return path;

            int startTri = navMesh->findTriangleContaining(start);
            int goalTri = navMesh->findTriangleContaining(goal);

            if (startTri == -1 || goalTri == -1) {
                return path; 
            }

            if (startTri == goalTri) {
                path.push_back(startTri);
                return path;
            }

            std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> openSet;
            std::unordered_map<int, double> gScores;
            std::unordered_map<int, int> cameFrom;

            double startH = heuristic(
                navMesh->getTriangle(startTri).tri.centroid(),
                navMesh->getTriangle(goalTri).tri.centroid()
            );

            openSet.push(AStarNode(startTri, 0, startH, -1));
            gScores[startTri] = 0;

            while (!openSet.empty()) {
                AStarNode current = openSet.top();
                openSet.pop();

                if (current.gCost > gScores[current.triangleId]) continue;

                if (current.triangleId == goalTri) {
                    int curr = goalTri;
                    while (curr != -1) {
                        path.push_back(curr);
                        curr = (cameFrom.find(curr) != cameFrom.end()) ? cameFrom[curr] : -1;
                    }
                    std::reverse(path.begin(), path.end());
                    return path;
                }

                const NavTriangle<T>& currentTri = navMesh->getTriangle(current.triangleId);

                for (int neighborId : currentTri.neighbors) {
                    Point<T> currentCenter = currentTri.tri.centroid();
                    Point<T> neighborCenter = navMesh->getTriangle(neighborId).tri.centroid();

                    double tentativeG = current.gCost + heuristic(currentCenter, neighborCenter);

                    if (gScores.find(neighborId) == gScores.end() ||
                        tentativeG < gScores[neighborId]) {

                        gScores[neighborId] = tentativeG;
                        cameFrom[neighborId] = current.triangleId;

                        double h = heuristic(
                            neighborCenter,
                            navMesh->getTriangle(goalTri).tri.centroid()
                        );

                        openSet.push(AStarNode(neighborId, tentativeG, h, current.triangleId));
                    }
                }
            }

            return path;
        }
    };

} // namespace Geom