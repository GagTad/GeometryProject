#pragma once

#include <optional>
#include "vector.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "polygon.hpp"
#include <vector>    
#include <algorithm>  
#include "ray.hpp"

namespace geom {

    template<size_t Dim, typename T>
    T distance(const Point<Dim, T>& p, const Line<Dim, T>& l) {
        Point<Dim, T> projected_point = l.project(p);

        return (p - projected_point).length();
    }

    template<typename T>
    T side_of_line(const Point<2, T>& p, const Line<2, T>& l) {
        Point<2, T> p2 = l.origin() + l.direction();
        return (p.coords[0] - l.origin().coords[0]) * (p2.coords[1] - l.origin().coords[1]) -
            (p.coords[1] - l.origin().coords[1]) * (p2.coords[0] - l.origin().coords[0]);
    }

    template<typename T>
    struct IntersectionResult2D {
        enum class Status { INTERSECTING, PARALLEL, COINCIDENT };

        Status status;
        std::optional<Point<2, T>> point;
    };

    template<typename T>
    T cross_product(const Vector<2, T>& a, const Vector<2, T>& b) {
        return a[0] * b[1] - a[1] * b[0];
    }

    template<typename T>
    IntersectionResult2D<T> intersection(const Line<2, T>& l1, const Line<2, T>& l2) {
        using Result = IntersectionResult2D<T>;

        const auto& p1 = l1.origin();
        const auto& dir1 = l1.direction();
        const auto& p2 = l2.origin();
        const auto& dir2 = l2.direction();

        const T dir_cross = cross_product(dir1, dir2);
        const Vector<2, T> p_diff = p2 - p1;

        if (std::abs(dir_cross) < Coord<T>::Epsilon) {
            if (std::abs(cross_product(p_diff, dir1)) < Coord<T>::Epsilon) {
                return { Result::Status::COINCIDENT, std::nullopt };
            }
            else {
                return { Result::Status::PARALLEL, std::nullopt };
            }
        }
        const T t = cross_product(p_diff, dir2) / dir_cross;
        const Point<2, T> intersection_point = p1 + dir1 * t;
        return { Result::Status::INTERSECTING, intersection_point };
    }

    template<typename T>
    struct SegmentIntersectionResult2D {
        enum class Status { NO_INTERSECTION, INTERSECTING, OVERLAPPING };

        Status status;
        std::optional<Point<2, T>> point;   
        std::optional<Segment<2, T>> segment;
    };

    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Segment<2, T>& s1, const Segment<2, T>& s2) {
        using Result = SegmentIntersectionResult2D<T>;

        Line<2, T> l1 = Line<2, T>::from_points(s1.p1(), s1.p2());
        Line<2, T> l2 = Line<2, T>::from_points(s2.p1(), s2.p2());

        auto line_isect = intersection(l1, l2);

        switch (line_isect.status) {
        case IntersectionResult2D<T>::Status::PARALLEL: {
            return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt };
        }

        case IntersectionResult2D<T>::Status::INTERSECTING: {
            Point<2, T> p = line_isect.point.value();
            if (contains(p, s1) && contains(p, s2)) {
                return { Result::Status::INTERSECTING, p, std::nullopt };
            }
            else {
                return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT: {
            auto p1 = s1.p1(), p2 = s1.p2();
            auto q1 = s2.p1(), q2 = s2.p2();
            if (p1[0] > p2[0]) std::swap(p1, p2);
            if (q1[0] > q2[0]) std::swap(q1, q2);

            if (p1[0] <= q2[0] && q1[0] <= p2[0]) {
                Point<2, T> overlap_start = (p1[0] > q1[0]) ? p1 : q1;
                Point<2, T> overlap_end = (p2[0] < q2[0]) ? p2 : q2;

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
        return { Result::Status::NO_INTERSECTION, std::nullopt, std::nullopt }; 
    }

    template<typename T>
    std::optional<Polygon<2, T>> convex_hull(std::vector<Point<2, T>>& points) {
        if (points.size() < 3) {
            return std::nullopt; 
        }

        std::sort(points.begin(), points.end(), [](const Point<2, T>& a, const Point<2, T>& b) {
            if (a[0] != b[0]) return a[0] < b[0];
            return a[1] < b[1];
            });

        std::vector<Point<2, T>> lower_hull;
        std::vector<Point<2, T>> upper_hull;

        for (const auto& p : points) {
            while (lower_hull.size() >= 2 &&
                cross_product(lower_hull.back() - lower_hull[lower_hull.size() - 2], p - lower_hull.back()) <= 0) {
                lower_hull.pop_back();
            }
            lower_hull.push_back(p);
        }

        for (int i = points.size() - 1; i >= 0; --i) {
            const auto& p = points[i];
            while (upper_hull.size() >= 2 &&
                cross_product(upper_hull.back() - upper_hull[upper_hull.size() - 2], p - upper_hull.back()) <= 0) {
                upper_hull.pop_back();
            }
            upper_hull.push_back(p);
        }

        upper_hull.pop_back();
        upper_hull.erase(upper_hull.begin());

        lower_hull.insert(lower_hull.end(), upper_hull.begin(), upper_hull.end());

        return Polygon<2, T>(lower_hull);
    }

    template<typename T>
    bool contains(const Point<2, T>& p, const Polygon<2, T>& polygon) {
        bool is_inside = false;
        const size_t num_verts = polygon.num_vertices();

        for (size_t i = 0; i < num_verts; ++i) {
            const auto& p1 = polygon.vertices()[i];
            const auto& p2 = polygon.vertices()[(i + 1) % num_verts];

            if (contains(p, Segment<2, T>(p1, p2))) {
                return true;
            }

            if (((p1[1] > p[1]) != (p2[1] > p[1]))) {
                const T x_intersection = (p2[0] - p1[0]) * (p[1] - p1[1]) / (p2[1] - p1[1]) + p1[0];

                if (x_intersection > p[0]) {
                    is_inside = !is_inside; // Փոխում ենք վիճակը (զույգ -> կենտ, կենտ -> զույգ)
                }
            }
        }

        return is_inside;
    }


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
            return { Result::Status::OVERLAPPING, std::nullopt, segment };
        }
        return { Result::Status::NO_INTERSECTION };
    }

    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Segment<2, T>& segment, const Line<2, T>& line) {
        return intersection(line, segment);
    }

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
            if (ray.contains(p)) { 
                return { Result::Status::INTERSECTING, p };
            }
            else {
                return { Result::Status::NO_INTERSECTION };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT:
            return { Result::Status::OVERLAPPING }; 
        }
        return { Result::Status::NO_INTERSECTION };
    }

    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Ray<2, T>& ray, const Line<2, T>& line) {
        return intersection(line, ray);
    }

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
            if (contains(p, segment) && ray.contains(p)) {
                return { Result::Status::INTERSECTING, p };
            }
            else {
                return { Result::Status::NO_INTERSECTION };
            }
        }

        case IntersectionResult2D<T>::Status::COINCIDENT: {
            auto p1 = segment.p1(), p2 = segment.p2();
            auto q1 = ray.origin();

            int axis = (std::abs(p2[0] - p1[0]) > std::abs(p2[1] - p1[1])) ? 0 : 1;
            if (p1[axis] > p2[axis]) std::swap(p1, p2);

            Point<2, T> overlap_start = (p1[axis] > q1[axis]) ? p1 : q1;
            Point<2, T> overlap_end = p2;

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

    template<typename T>
    SegmentIntersectionResult2D<T> intersection(const Ray<2, T>& ray, const Segment<2, T>& segment) {
        return intersection(segment, ray);
    }

} // namespace geom

