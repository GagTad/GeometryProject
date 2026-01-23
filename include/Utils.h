#pragma once
#include <cmath>
#include <type_traits>
#include <algorithm>

namespace Geom {

    const double EPS = 1e-9;

    inline bool equals(double a, double b) {
        return std::abs(a - b) < EPS;
    }

    inline int sign(double x) {
        if (x < -EPS) return -1;
        if (x > EPS) return 1;
        return 0;
    }
}
