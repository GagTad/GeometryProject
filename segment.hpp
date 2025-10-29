#pragma once

#pragma once

#include "vector.hpp"
#include <algorithm> // std::clamp-ի համար (C++17)

namespace geom {

    template<size_t Dim, typename T>
    class Segment {
    public:
        using point_type = Point<Dim, T>;
		using vector_type = Vector<Dim, T>;
    private:
        point_type m_p1;
        point_type m_p2;

    public:
        // --- Կոնստրուկտոր ---
        Segment(const point_type& p1, const point_type& p2) : m_p1(p1), m_p2(p2) {
            // Կարող ենք ստուգում ավելացնել, որ կետերը չհամընկնեն,
            // բայց զրոյական երկարությամբ հատվածը նույնպես թույլատրելի է։
        }

        // --- Հասանելիության մեթոդներ (Getters) ---
        const point_type& p1() const { return m_p1; }
        const point_type& p2() const { return m_p2; }

        // --- Երկրաչափական մեթոդներ ---

        // Հատվածի երկարության քառակուսին
        T length_sq() const {
            return (m_p2 - m_p1).length_sq();
        }

        // Հատվածի երկարությունը
        T length() const {
            return (m_p2 - m_p1).length();
        }

        // Գտնում է p կետի պրոյեկցիան հատվածի վրա
        point_type project(const point_type& p) const {
            const vector_type v = m_p2 - m_p1;
            const T len_sq = v.length_sq();

            // Եթե հատվածը կետ է (զրո երկարություն)
            if (len_sq < Coord<T>::Epsilon) {
                return m_p1;
            }

            // t-ն ցույց է տալիս, թե p-ի պրոյեկցիան որքան հեռու է p1-ից p2 ուղղությամբ
            // t = dot_product(AP, AB) / |AB|^2
            const T t = dot_product(p - m_p1, v) / len_sq;

            // Եթե t < 0, պրոյեկցիան p1-ից "ձախ" է, ամենամոտ կետը p1-ն է։
            // Եթե t > 1, պրոյեկցիան p2-ից "աջ" է, ամենամոտ կետը p2-ն է։
            // Եթե 0 <= t <= 1, պրոյեկցիան հատվածի վրա է։
            const T clamped_t = std::clamp(t, T(0), T(1));

            return m_p1 + v * clamped_t;
        }
    };

    // --- Հարակից ֆունկցիաներ (դրսում) ---

    // Հաշվում է p կետի հեռավորությունը s հատվածից
    template<size_t Dim, typename T>
    T distance(const Point<Dim, T>& p, const Segment<Dim, T>& s) {
        const Point<Dim, T> projected_point = s.project(p);
        return (p - projected_point).length();
    }

    // Ստուգում է՝ արդյոք p կետը գտնվում է s հատվածի վրա
    template<size_t Dim, typename T>
    bool contains(const Point<Dim, T>& p, const Segment<Dim, T>& s) {
        // Կետը հատվածի վրա է, եթե նրա հեռավորությունը հատվածից զրո է (epsilon-ի ճշտությամբ)
        return distance(p, s) < Coord<T>::Epsilon;
    }


    // --- Տիպերի այլանուններ ---
    template<typename T> using Segment2 = Segment<2, T>;
    template<typename T> using Segment3 = Segment<3, T>;

    using Segment2d = Segment2<double>;
    using Segment3d = Segment3<double>;

} // namespace geom