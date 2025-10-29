#pragma once

#include <cmath>
#include <type_traits>
#include <iostream>

namespace geom {

    template<typename T>
    struct Coord {
        T value;

        static_assert(std::is_floating_point<T>::value, "Coord<T> is intended for floating-point types like double or float.");

        // Համեմատության համար օգտագործվող թույլատրելի սխալանք (epsilon)
        static constexpr T Epsilon = 1e-9;

        // --- Կոնստրուկտորներ ---
        Coord() : value(0) {}
        Coord(const T& val) : value(val) {} // թույլ է տալիս գրել Coord c = 5.0;

        // --- Համեմատության օպերատորներ ---
        bool operator==(const Coord& other) const { return std::abs(value - other.value) <= Epsilon; }
        bool operator!=(const Coord& other) const { return !(*this == other); }
        bool operator<(const Coord& other) const { return value < other.value && !(*this == other); }
        bool operator>(const Coord& other) const { return value > other.value && !(*this == other); }
        bool operator<=(const Coord& other) const { return value < other.value || (*this == other); }
        bool operator>=(const Coord& other) const { return value > other.value || (*this == other); }

        // --- Փոխակերպում դեպի T ---
        // թույլ է տալիս օգտագործել Coord-ը այնտեղ, որտեղ T է սպասվում (օր.՝ std::sqrt(my_coord))
        operator T() const { return value; }

        // --- Միանիշ (unary) օպերատորներ ---
        Coord operator-() const { return Coord(-value); }

        // --- Թվաբանական վերագրման օպերատորներ ---
        Coord& operator+=(const Coord& other) { value += other.value; return *this; }
        Coord& operator-=(const Coord& other) { value -= other.value; return *this; }
        Coord& operator*=(const Coord& other) { value *= other.value; return *this; }
        Coord& operator/=(const Coord& other) { value /= other.value; return *this; }
    };

    // --- Երկնիշ (binary) թվաբանական օպերատորներ ---
    // Սրանք թույլ են տալիս գրել c3 = c1 + c2;
    template<typename T> Coord<T> operator+(Coord<T> lhs, const Coord<T>& rhs) { lhs += rhs; return lhs; }
    template<typename T> Coord<T> operator-(Coord<T> lhs, const Coord<T>& rhs) { lhs -= rhs; return lhs; }
    template<typename T> Coord<T> operator*(Coord<T> lhs, const Coord<T>& rhs) { lhs *= rhs; return lhs; }
    template<typename T> Coord<T> operator/(Coord<T> lhs, const Coord<T>& rhs) { lhs /= rhs; return lhs; }

    // --- Տպելու համար օպերատոր ---
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const Coord<T>& c) {
        os << c.value;
        return os;
    }

} // namespace geom