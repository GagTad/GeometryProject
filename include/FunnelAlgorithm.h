#pragma once
#include "NavMesh.h"
#include <vector>

namespace Geom {

    template<typename T>
    class FunnelAlgorithm {
    private:
        const NavMesh<T>* navMesh;

        T triArea2(const Point<T>& a, const Point<T>& b, const Point<T>& c) const {
            return crossProduct(b - a, c - a);
        }

    public:
        FunnelAlgorithm(const NavMesh<T>* mesh) : navMesh(mesh) {}

        std::vector<Point<T>> smoothPath(const std::vector<int>& trianglePath,
            const Point<T>& start,
            const Point<T>& goal) {
            std::vector<Point<T>> smoothed;

            if (trianglePath.empty()) return smoothed;

            if (trianglePath.size() == 1) {
                smoothed.push_back(start);
                smoothed.push_back(goal);
                return smoothed;
            }

            struct Portal {
                Point<T> left, right;
            };

            std::vector<Portal> portals;

            for (size_t i = 0; i < trianglePath.size() - 1; i++) {
                Point<T> edgeStart, edgeEnd;
                if (navMesh->getSharedEdge(trianglePath[i], trianglePath[i + 1],
                    edgeStart, edgeEnd)) {
                    Portal p;

                    const Triangle<T>& tri = navMesh->getTriangle(trianglePath[i]).tri;
                    Point<T> triCenter = tri.centroid();

                    if (triArea2(edgeStart, edgeEnd, triCenter) > 0) {
                        p.left = edgeEnd;
                        p.right = edgeStart;
                    }
                    else {
                        p.left = edgeStart;
                        p.right = edgeEnd;
                    }

                    portals.push_back(p);
                }
            }

            Portal finalPortal;
            finalPortal.left = goal;
            finalPortal.right = goal;
            portals.push_back(finalPortal);

            if (portals.empty()) {
                smoothed.push_back(start);
                smoothed.push_back(goal);
                return smoothed;
            }

            Point<T> apex = start;
            Point<T> leftPoint = portals[0].left;
            Point<T> rightPoint = portals[0].right;
            int apexIndex = 0;
            int leftIndex = 0;
            int rightIndex = 0;

            smoothed.push_back(apex);

            for (size_t i = 0; i < portals.size(); i++) {
                Point<T> newLeft = portals[i].left;
                Point<T> newRight = portals[i].right;


                if (triArea2(apex, rightPoint, newRight) >= -EPS) {

                    if (apex == rightPoint || triArea2(apex, leftPoint, newRight) <= EPS) {
                        rightPoint = newRight;
                        rightIndex = i;
                    }
                    else {
                        apex = leftPoint;
                        smoothed.push_back(apex);
                        apexIndex = leftIndex;
                        leftPoint = apex;
                        rightPoint = apex;
                        leftIndex = apexIndex;
                        rightIndex = apexIndex;
                        i = apexIndex;
                        continue;
                    }
                }


                if (triArea2(apex, leftPoint, newLeft) <= EPS) {

                    if (apex == leftPoint || triArea2(apex, rightPoint, newLeft) >= -EPS) {
                        leftPoint = newLeft;
                        leftIndex = i;
                    }
                    else {
                        apex = rightPoint;
                        smoothed.push_back(apex);

                        apexIndex = rightIndex;
                        leftPoint = apex;
                        rightPoint = apex;
                        leftIndex = apexIndex;
                        rightIndex = apexIndex;

                        i = apexIndex;
                        continue;
                    }
                }
            }

            if (!(smoothed.back() == goal)) {
                smoothed.push_back(goal);
            }

            return smoothed;
        }
    };

} // namespace Geom