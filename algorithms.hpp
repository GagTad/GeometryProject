#pragma once

#include <optional>
#include "vector.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "polygon.hpp"
#include <vector>     // std::vector-ի համար
#include <algorithm>  // std::sort-ի համար
#include "ray.hpp"
namespace geom {

    // --- Point-Line Algorithms ---

    // 1. Հաշվում է p կետի հեռավորությունը l ուղղից
    template<size_t Dim, typename T>
    T distance(const Point<Dim, T>& p, const Line<Dim, T>& l) {
        // Գտնում ենք կետի պրոյեկցիան ուղղի վրա
        Point<Dim, T> projected_point = l.project(p);

        // Հեռավորությունը կետի և իր պրոյեկցիայի միջև եղած հեռավորությունն է
        return (p - projected_point).length();
    }

    // 2. (Միայն 2D-ի համար) Որոշում է, թե p կետը ուղղի որ կողմում է գտնվում
    // Վերադարձնում է.
    // > 0 եթե մի կողմում է
    // < 0 եթե մյուս կողմում է
    // = 0 (epsilon-ի սահմաններում) եթե ուղղի վրա է
    template<typename T>
    T side_of_line(const Point<2, T>& p, const Line<2, T>& l) {
        // Օգտագործում ենք վեկտորական արտադրյալի 2D անալոգը (cross product z-component)
        // l.origin()-ից p-ին և l.origin()-ից l.direction()-ին տարած վեկտորների համար
        Point<2, T> p2 = l.origin() + l.direction();
        return (p.coords[0] - l.origin().coords[0]) * (p2.coords[1] - l.origin().coords[1]) -
            (p.coords[1] - l.origin().coords[1]) * (p2.coords[0] - l.origin().coords[0]);
    }

    // 2D ուղիղների հատման արդյունքը նկարագրող կառուցվածք
    template<typename T>
    struct IntersectionResult2D {
        enum class Status { INTERSECTING, PARALLEL, COINCIDENT };

        Status status;
        std::optional<Point<2, T>> point;
    };

    // 2D վեկտորների վեկտորական արտադրյալի անալոգ (cross product)
    // Վերադարձնում է Z կոմպոնենտի մեծությունը
    template<typename T>
    T cross_product(const Vector<2, T>& a, const Vector<2, T>& b) {
        // Ուղղված է՝ օգտագործում ենք public operator[]-ը coords-ի փոխարեն
        return a[0] * b[1] - a[1] * b[0];
    }

    // Հաշվում է երկու 2D ուղիղների հատումը
    template<typename T>
    IntersectionResult2D<T> intersection(const Line<2, T>& l1, const Line<2, T>& l2) {
        using Result = IntersectionResult2D<T>;

        const auto& p1 = l1.origin();
        const auto& dir1 = l1.direction();
        const auto& p2 = l2.origin();
        const auto& dir2 = l2.direction();

        const T dir_cross = cross_product(dir1, dir2);
        const Vector<2, T> p_diff = p2 - p1;

        // Ստուգում ենք՝ արդյոք ուղղորդող վեկտորները զուգահեռ են
        if (std::abs(dir_cross) < Coord<T>::Epsilon) {
            // Ուղիղները զուգահեռ են։ Ստուգենք՝ համընկնո՞ւմ են, թե՞ ոչ։
            // Եթե p1-ից p2 տարած վեկտորը նույնպես զուգահեռ է, ուրեմն համընկնում են։
            if (std::abs(cross_product(p_diff, dir1)) < Coord<T>::Epsilon) {
                return { Result::Status::COINCIDENT, std::nullopt };
            }
            else {
                return { Result::Status::PARALLEL, std::nullopt };
            }
        }

        // Ուղիղները հատվում են։ Գտնենք t պարամետրը l1-ի համար։
        // t = (p2 - p1) x dir2 / (dir1 x dir2)
        const T t = cross_product(p_diff, dir2) / dir_cross;

        const Point<2, T> intersection_point = p1 + dir1 * t;

        return { Result::Status::INTERSECTING, intersection_point };
    }

    // ... intersection(Line, Line) ֆունկցիայից հետո ...

// --- Segment-Segment Algorithms ---

// 2D հատվածների հատման արդյունքը նկարագրող կառուցվածք
    template<typename T>
    struct SegmentIntersectionResult2D {
        enum class Status { NO_INTERSECTION, INTERSECTING, OVERLAPPING };

        Status status;
        std::optional<Point<2, T>> point;      // INTERSECTING դեպքի համար
        std::optional<Segment<2, T>> segment; // OVERLAPPING դեպքի համար
    };

    // Հաշվում է երկու 2D հատվածների հատումը
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Segment<2, T>& s1, const Segment<2, T>& s2) {
        using Result = SegmentIntersectionResult2D<T>;

        // 1. "Լայն" փուլ. ստուգում ենք պարունակող ուղիղները
        Line<2, T> l1 = Line<2, T>::from_points(s1.p1(), s1.p2());
        Line<2, T> l2 = Line<2, T>::from_points(s2.p1(), s2.p2());

        auto line_isect = intersection(l1, l2);

        switch (line_isect.status) {
        case IntersectionResult2D<T>::Status::PARALLEL: {
            // Եթե ուղիղները զուգահեռ են, հատվածները չեն հատվում
            return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt };
        }

        case IntersectionResult2D<T>::Status::INTERSECTING: {
            // 2. "Նեղ" փուլ. ստուգում ենք՝ արդյոք հատման կետը երկու հատվածների վրա է
            Point<2, T> p = line_isect.point.value();
            if (contains(p, s1) && contains(p, s2)) {
                return { Result::Status::INTERSECTING, p, std::nullopt };
            }
            else {
                return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT: {
            // 3. Ամենաբարդ դեպքը. հատվածները նույն ուղղի վրա են
            // Ստուգում ենք 1D վերադրումը (overlap)
            auto p1 = s1.p1(), p2 = s1.p2();
            auto q1 = s2.p1(), q2 = s2.p2();

            // Դասավորում ենք յուրաքանչյուր հատվածի կետերը (օրինակ՝ ըստ x-ի)
            if (p1[0] > p2[0]) std::swap(p1, p2);
            if (q1[0] > q2[0]) std::swap(q1, q2);

            // Ստուգում ենք, որ երկու [p1, p2] և [q1, q2] 1D հատվածները վերադրվում են
            if (p1[0] <= q2[0] && q1[0] <= p2[0]) {
                // Վերադրում կա։ Գտնենք վերադրման հատվածի ծայրակետերը։
                Point<2, T> overlap_start = (p1[0] > q1[0]) ? p1 : q1;
                Point<2, T> overlap_end = (p2[0] < q2[0]) ? p2 : q2;

                // Եթե ծայրակետերը համարյա նույնն են, ուրեմն մեկ կետում են հպվում
                if ((overlap_start - overlap_end).length_sq() < Coord<T>::Epsilon) {
                    return { Result::Status::INTERSECTING, overlap_start, std::nullopt };
                }

                return { Result::Status::OVERLAPPING, std::nullopt, Segment<2, T>(overlap_start, overlap_end) };
            }
            else {
                return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt };
            }
        }
        }
        return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt }; // Just in case
    }

    // Գտնում է 2D կետերի բազմության ուռուցիկ բազմությունը (Convex Hull)
    template<typename T>
    std::optional<Polygon<2, T>> convex_hull(std::vector<Point<2, T>>& points) {
        if (points.size() < 3) {
            return std::nullopt; // Ուռուցիկ բազմությունը իմաստ ունի առնվազն 3 կետի համար
        }

        // 1. Սորտավորում
        std::sort(points.begin(), points.end(), [](const Point<2, T>& a, const Point<2, T>& b) {
            if (a[0] != b[0]) return a[0] < b[0];
            return a[1] < b[1];
            });

        std::vector<Point<2, T>> lower_hull;
        std::vector<Point<2, T>> upper_hull;

        // 2. Ստորին կեղևի կառուցում
        for (const auto& p : points) {
            while (lower_hull.size() >= 2 &&
                cross_product(lower_hull.back() - lower_hull[lower_hull.size() - 2], p - lower_hull.back()) <= 0) {
                lower_hull.pop_back();
            }
            lower_hull.push_back(p);
        }

        // 3. Վերին կեղևի կառուցում (հակառակ ուղղությամբ)
        for (int i = points.size() - 1; i >= 0; --i) {
            const auto& p = points[i];
            while (upper_hull.size() >= 2 &&
                cross_product(upper_hull.back() - upper_hull[upper_hull.size() - 2], p - upper_hull.back()) <= 0) {
                upper_hull.pop_back();
            }
            upper_hull.push_back(p);
        }

        // 4. Միավորում
        // Հեռացնում ենք վերին կեղևի առաջին և վերջին կետերը, քանի որ դրանք
        // կրկնում են ստորին կեղևի վերջին և առաջին կետերին
        upper_hull.pop_back();
        upper_hull.erase(upper_hull.begin());

        lower_hull.insert(lower_hull.end(), upper_hull.begin(), upper_hull.end());

        return Polygon<2, T>(lower_hull);
    }
    // ... convex_hull ֆունկցիայից հետո ...

// --- Point-Polygon Algorithms ---

// Ստուգում է, թե արդյոք p կետը գտնվում է polygon-ի ներսում (Ray Casting ալգորիթմ)
    template<typename T>
    bool contains(const Point<2, T>& p, const Polygon<2, T>& polygon) {
        bool is_inside = false;
        const size_t num_verts = polygon.num_vertices();

        for (size_t i = 0; i < num_verts; ++i) {
            const auto& p1 = polygon.vertices()[i];
            const auto& p2 = polygon.vertices()[(i + 1) % num_verts];

            // Նախնական ստուգում. եթե կետը եզրագծի վրա է, համարում ենք, որ ներսում է
            if (contains(p, Segment<2, T>(p1, p2))) {
                return true;
            }

            // Ray Casting տրամաբանությունը
            // Ստուգում ենք, որ կողմի y կոորդինատները գտնվում են ճառագայթի տարբեր կողմերում
            if (((p1[1] > p[1]) != (p2[1] > p[1]))) {
                // Հաշվում ենք, թե հորիզոնական ճառագայթը որ x կոորդինատում է հատում կողմը
                const T x_intersection = (p2[0] - p1[0]) * (p[1] - p1[1]) / (p2[1] - p1[1]) + p1[0];

                // Եթե հատման կետը գտնվում է մեր կետից աջ, ուրեմն ճառագայթը հատել է կողմը
                if (x_intersection > p[0]) {
                    is_inside = !is_inside; // Փոխում ենք վիճակը (զույգ -> կենտ, կենտ -> զույգ)
                }
            }
        }

        return is_inside;
    }
    // ... intersection(Segment, Segment) ֆունկցիայից հետո ...

// --- Mixed Intersection Algorithms ---

// Ուղղի և հատվածի հատումը
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Line<2, T>& line, const Segment<2, T>& segment) {
        using Result = SegmentIntersectionResult2D<T>;

        Line<2, T> seg_line = Line<2, T>::from_points(segment.p1(), segment.p2());
        auto line_isect = intersection(line, seg_line);

        switch (line_isect.status) {
        case IntersectionResult2D<T>::Status::PARALLEL:
            return { Result::Status::NO_INTERSECTION };

        case IntersectionResult2D<T>::Status::INTERSECTING: {
            Point<2, T> p = line_isect.point.value();
            if (contains(p, segment)) {
                return { Result::Status::INTERSECTING, p };
            }
            else {
                return { Result::Status::NO_INTERSECTION };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT:
            // Եթե ուղիղները համընկնում են, ապա հատումը հենց հատվածն է
            return { Result::Status::OVERLAPPING, std::nullopt, segment };
        }
        return { Result::Status::NO_INTERSECTION };
    }

    // Հարմարության համար՝ նույն ֆունկցիան հակառակ արգումենտներով
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Segment<2, T>& segment, const Line<2, T>& line) {
        return intersection(line, segment);
    }

    // ... intersection(Segment, Line) ֆունկցիայից հետո ...

// Ուղղի և ճառագայթի հատումը
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Line<2, T>& line, const Ray<2, T>& ray) {
        using Result = SegmentIntersectionResult2D<T>;

        Line<2, T> ray_line = Line<2, T>::from_point_direction(ray.origin(), ray.direction());
        auto line_isect = intersection(line, ray_line);

        switch (line_isect.status) {
        case IntersectionResult2D<T>::Status::PARALLEL:
            return { Result::Status::NO_INTERSECTION };

        case IntersectionResult2D<T>::Status::INTERSECTING: {
            Point<2, T> p = line_isect.point.value();
            if (ray.contains(p)) { // Ստուգում ենք միայն ճառագայթի համար
                return { Result::Status::INTERSECTING, p };
            }
            else {
                return { Result::Status::NO_INTERSECTION };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT:
            // Եթե ուղիղները համընկնում են, հատումը հենց ճառագայթն է
            // Մեր Result struct-ը չի կարող Ray պահել, բայց կարող ենք նշել որպես OVERLAPPING
            return { Result::Status::OVERLAPPING }; // Նշում է, որ հատումը անվերջ է
        }
        return { Result::Status::NO_INTERSECTION };
    }

    // Հակառակ արգումենտներով տարբերակը
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Ray<2, T>& ray, const Line<2, T>& line) {
        return intersection(line, ray);
    }
    // (Այստեղ հետագայում կավելացնենք մյուս զույգերի համար ֆունկցիաները)
    // ... intersection(Ray, Line) ֆունկցիայից հետո ...

// Հատվածի և ճառագայթի հատումը
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Segment<2, T>& segment, const Ray<2, T>& ray) {
        using Result = SegmentIntersectionResult2D<T>;

        Line<2, T> seg_line = Line<2, T>::from_points(segment.p1(), segment.p2());
        Line<2, T> ray_line = Line<2, T>::from_point_direction(ray.origin(), ray.direction());
        auto line_isect = intersection(seg_line, ray_line);

        switch (line_isect.status) {
        case IntersectionResult2D<T>::Status::PARALLEL:
            return { Result::Status::NO_INTERSECTION };

        case IntersectionResult2D<T>::Status::INTERSECTING: {
            Point<2, T> p = line_isect.point.value();
            // Ստուգում ենք, որ կետը գտնվում է ԵՎ՛ հատվածի, ԵՎ՛ ճառագայթի վրա
            if (contains(p, segment) && ray.contains(p)) {
                return { Result::Status::INTERSECTING, p };
            }
            else {
                return { Result::Status::NO_INTERSECTION };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT: {
            // Ստուգում ենք 1D վերադրումը հատվածի [p1,p2] և ճառագայթի [q1, +inf) միջև
            auto p1 = segment.p1(), p2 = segment.p2();
            auto q1 = ray.origin();

            int axis = (std::abs(p2[0] - p1[0]) > std::abs(p2[1] - p1[1])) ? 0 : 1;
            if (p1[axis] > p2[axis]) std::swap(p1, p2);

            // Գտնում ենք վերադրման սկիզբը և վերջը
            Point<2, T> overlap_start = (p1[axis] > q1[axis]) ? p1 : q1;
            Point<2, T> overlap_end = p2;

            // Ստուգում ենք՝ արդյոք ճառագայթը ճիշտ ուղղությամբ է շարժվում
            bool ray_goes_towards_segment = dot_product(p2 - q1, ray.direction()) >= 0;

            if (ray_goes_towards_segment && overlap_start[axis] <= overlap_end[axis]) {
                if ((overlap_start - overlap_end).length_sq() < Coord<T>::Epsilon * Coord<T>::Epsilon) {
                    return { Result::Status::INTERSECTING, overlap_start };
                }
                return { Result::Status::OVERLAPPING, std::nullopt, Segment<2, T>(overlap_start, overlap_end) };
            }
            return { Result::Status::NO_INTERSECTION };
        }
        }
        return { Result::Status::NO_INTERSECTION };
    }

    // Հակառակ արգումենտներով տարբերակը
    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Ray<2, T>& ray, const Segment<2, T>& segment) {
        return intersection(segment, ray);
    }
} // namespace geom