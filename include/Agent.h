#pragma once
#include "geometry.h"
#include <vector>

namespace Geom {

    class Agent {
    private:
        Point<double> position;
        std::vector<Point<double>> path;
        int currentTargetIndex;
        double speed;
        bool active;

    public:
        Agent() : speed(150.0), currentTargetIndex(0), active(false) {}

        void setPosition(const Point<double>& p) {
            position = p;
        }

        Point<double> getPosition() const {
            return position;
        }

        void setPath(const std::vector<Point<double>>& newPath) {
            path = newPath;
            if (path.size() > 0) {
                position = path[0];
                currentTargetIndex = 1;
                active = true;
            }
        }

        void stop() {
            active = false;
        }

        void update(double dt) {
            if (!active || path.empty() || currentTargetIndex >= path.size()) {
                return;
            }

            Point<double> target = path[currentTargetIndex];
            Vector<double> direction = target - position;
            double distToTarget = length(direction);

            if (distToTarget < 5.0) { 
                currentTargetIndex++;
                if (currentTargetIndex >= path.size()) {
                    active = false;
                }
                return;
            }

            Vector<double> dirNorm = normalize(direction);

            position.x += dirNorm.dx * speed * dt;
            position.y += dirNorm.dy * speed * dt;
        }

        bool isActive() const { return active; }
    };
}