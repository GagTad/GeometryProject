#include <iostream>
#include "vector.hpp"
#include "line.hpp"   // *** ՆՈՐ ԸՆԴԳՐԿՈՒՄ ***
#include "algorithms.hpp"
#include "segment.hpp"
#include "polygon.hpp"
#include "ray.hpp"

// Նախորդ թեստերը տեղափոխում ենք իրենց ֆունկցիայի մեջ
void run_vector_tests() {
    std::cout << "--- Running Vector/Point Tests ---" << std::endl;

    geom::Point2d p1(3.0, 4.0);
    geom::Vector2d v1(1.0, -2.0);
    std::cout << "Test 1.1: Point p1 = " << p1 << std::endl;
    std::cout << "Test 1.2: Vector v1 = " << v1 << std::endl;

    geom::Point2d p2 = p1 + v1;
    std::cout << "Test 2.1: p1 + v1 = " << p2 << std::endl;

    geom::Vector2d v2 = p1 - geom::Point2d(1.0, 1.0);
    std::cout << "Test 2.2: p1 - (1,1) = " << v2 << std::endl;

    geom::Vector2d v3 = v1 * 2.5;
    std::cout << "Test 3.1: v1 * 2.5 = " << v3 << std::endl;

    std::cout << "Test 4.1: Length of p1 = " << p1.length() << std::endl;

    geom::Vector2d unit_x(1, 0), unit_y(0, 1);
    std::cout << "Test 5.1: dot_product((1,0), (0,1)) = " << dot_product(unit_x, unit_y) << std::endl;
    std::cout << "Test 5.2: dot_product(p1, p1) = " << dot_product(p1, p1) << " (should be 25)" << std::endl;

    geom::Point2d p_a(1.0, 2.0), p_b(1.0000000001, 2.0), p_c(1.001, 2.0);

    std::cout << "Test 6.1: Comparing (1, 2) and (1.0000000001, 2)... ";
    if (p_a == p_b) std::cout << "SUCCESS (Equal)" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 6.2: Comparing (1, 2) and (1.001, 2)... ";
    if (p_a != p_c) std::cout << "SUCCESS (Not Equal)" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "--- Vector/Point Tests Finished ---\n" << std::endl;
}

// *** ՆՈՐ ՖՈՒՆԿՑԻԱ LINE-Ի ՀԱՄԱՐ ***
void run_line_tests() {
    std::cout << "--- Running Line Tests ---" << std::endl;

    geom::Point2d p1(2, 2);
    geom::Point2d p2(6, 5); // Vector p2-p1 is (4, 3), length is 5
    geom::Line2d line = geom::Line2d::from_points(p1, p2);


    // Test 1: Construction and Normalization
    geom::Vector2d expected_dir(0.8, 0.6); // (4,3) / 5 = (0.8, 0.6)
    std::cout << "Test 1.1: Line direction normalization... ";
    if (line.direction() == expected_dir) {
        std::cout << "SUCCESS. Direction: " << line.direction() << std::endl;
    }
    else {
        std::cout << "FAILED. Expected " << expected_dir << " but got " << line.direction() << std::endl;
    }

    // Test 2: point_at(t)
    std::cout << "Test 2.1: point_at(0) should be the origin... ";
    if (line.point_at(0) == p1) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 2.2: point_at(5) should be p2... ";
    // t=5, because distance between p1 and p2 is 5, and direction is normalized.
    if (line.point_at(5) == p2) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED. Got " << line.point_at(5) << std::endl;

    // Test 3: project(p)
    geom::Point2d p_off_line(10, 0);
    geom::Point2d expected_proj(5, 3.5); // Calculated manually
    std::cout << "Test 3.1: Projecting point (10, 0)... ";
    if (line.project(p_off_line) == expected_proj) {
        std::cout << "SUCCESS. Projection: " << line.project(p_off_line) << std::endl;
    }
    else {
        std::cout << "FAILED. Expected " << expected_proj << " but got " << line.project(p_off_line) << std::endl;
    }

    // Test 4: contains(p)
    geom::Point2d p_on_line = line.point_at(-10); // A point definitely on the line
    geom::Point2d p_almost_on_line = p_on_line + geom::Vector2d(0, 1e-10); // Very close point
    geom::Point2d p_far_from_line(0, 0);

    std::cout << "Test 4.1: Contains point on line... ";
    if (line.contains(p_on_line)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 4.2: Contains point *almost* on line... ";
    if (line.contains(p_almost_on_line)) std::cout << "SUCCESS (Epsilon comparison works!)" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 4.3: Contains point far from line... ";
    if (!line.contains(p_far_from_line)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "--- Line Tests Finished ---" << std::endl;

    // Test 5: Point-Line distance
    std::cout << "Test 5.1: Distance from p1 (on line) to line... ";
    // p1 = (2,2), այն ուղղի վրա է, հեռավորությունը պետք է լինի 0
    if (geom::distance(p1, line) < geom::Coord<double>::Epsilon) {
        std::cout << "SUCCESS. Distance: " << geom::distance(p1, line) << std::endl;
    }
    else {
        std::cout << "FAILED. Distance: " << geom::distance(p1, line) << std::endl;
    }

    geom::Point2d p_for_dist(7, -2);
    // (7,-2) կետից մինչև (6.16, 5.12) պրոյեկցիան ունեցող կետ հեռավորությունը 7.2 է
    double expected_dist = 7.2;
    std::cout << "Test 5.2: Distance from (7, -2) to line... ";
    if (std::abs(geom::distance(p_for_dist, line) - expected_dist) < geom::Coord<double>::Epsilon) {
        std::cout << "SUCCESS. Distance: " << geom::distance(p_for_dist, line) << std::endl;
    }
    else {
        std::cout << "FAILED. Expected ~" << expected_dist << ", got " << geom::distance(p_for_dist, line) << std::endl;
    }
}


void run_segment_tests() {
    std::cout << "\n--- Running Segment Tests ---" << std::endl;

    geom::Point2d p1(0, 0);
    geom::Point2d p2(10, 0);
    geom::Segment2d seg(p1, p2);

    // Test 1: Length
    std::cout << "Test 1.1: Segment length... ";
    if (std::abs(seg.length() - 10.0) < geom::Coord<double>::Epsilon) {
        std::cout << "SUCCESS" << std::endl;
    }
    else {
        std::cout << "FAILED. Got " << seg.length() << std::endl;
    }

    // Test 2: Projection
    geom::Point2d p_inside(5, 5);
    geom::Point2d p_before(-5, 5);
    geom::Point2d p_after(15, 5);

    std::cout << "Test 2.1: Projecting internal point (5, 5)... ";
    if (seg.project(p_inside) == geom::Point2d(5, 0)) {
        std::cout << "SUCCESS. Got " << seg.project(p_inside) << std::endl;
    }
    else {
        std::cout << "FAILED" << std::endl;
    }

    std::cout << "Test 2.2: Projecting external point before p1... ";
    if (seg.project(p_before) == p1) {
        std::cout << "SUCCESS. Got " << seg.project(p_before) << std::endl;
    }
    else {
        std::cout << "FAILED" << std::endl;
    }

    std::cout << "Test 2.3: Projecting external point after p2... ";
    if (seg.project(p_after) == p2) {
        std::cout << "SUCCESS. Got " << seg.project(p_after) << std::endl;
    }
    else {
        std::cout << "FAILED" << std::endl;
    }

    // Test 3: Contains
    std::cout << "Test 3.1: Contains point in the middle... ";
    if (geom::contains(geom::Point2d(3, 0), seg)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 3.2: Does not contain point on supporting line... ";
    if (!geom::contains(geom::Point2d(12, 0), seg)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 3.3: Does not contain point off the line... ";
    if (!geom::contains(geom::Point2d(5, 5), seg)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "--- Segment Tests Finished ---" << std::endl;
}

void run_intersection_tests() {
    std::cout << "\n--- Running Intersection Tests ---" << std::endl;

    using Status = geom::IntersectionResult2D<double>::Status;

    // Case 1: Simple intersection
    geom::Line2d l1 = geom::Line2d::from_points(geom::Point2d(0, 1), geom::Point2d(10, 1)); // Horizontal line y=1
    geom::Line2d l2 = geom::Line2d::from_points(geom::Point2d(5, 0), geom::Point2d(5, 10)); // Vertical line x=5

    std::cout << "Test 1.1: Simple intersection... ";
    auto result1 = geom::intersection(l1, l2);
    if (result1.status == Status::INTERSECTING && result1.point.value() == geom::Point2d(5, 1)) {
        std::cout << "SUCCESS. Point: " << result1.point.value() << std::endl;
    }
    else {
        std::cout << "FAILED" << std::endl;
    }

    // Case 2: Parallel lines
    geom::Line2d l3 = geom::Line2d::from_points(geom::Point2d(0, 3), geom::Point2d(10, 3)); // Horizontal line y=3
    std::cout << "Test 1.2: Parallel lines... ";
    auto result2 = geom::intersection(l1, l3);
    if (result2.status == Status::PARALLEL) {
        std::cout << "SUCCESS" << std::endl;
    }
    else {
        std::cout << "FAILED" << std::endl;
    }

    // Case 3: Coincident lines
    geom::Line2d l4 = geom::Line2d::from_points(geom::Point2d(-5, 1), geom::Point2d(20, 1)); // Another line on y=1
    std::cout << "Test 1.3: Coincident lines... ";
    auto result3 = geom::intersection(l1, l4);
    if (result3.status == Status::COINCIDENT) {
        std::cout << "SUCCESS" << std::endl;
    }
    else {
        std::cout << "FAILED" << std::endl;
    }

    // Case 4: General intersection
    geom::Line2d l5 = geom::Line2d::from_points(geom::Point2d(2, 2), geom::Point2d(6, 5));   // y = 0.75x + 0.5
    geom::Line2d l6 = geom::Line2d::from_points(geom::Point2d(0, 5), geom::Point2d(6, 2));   // y = -0.5x + 5
    // Intersection point calculated manually: (3.6, 3.2)

    std::cout << "Test 1.4: General intersection... ";
    auto result4 = geom::intersection(l5, l6);
    if (result4.status == Status::INTERSECTING && result4.point.value() == geom::Point2d(3.6, 3.2)) {
        std::cout << "SUCCESS. Point: " << result4.point.value() << std::endl;
    }
    else {
        std::cout << "FAILED. ";
        if (result4.point.has_value()) {
            std::cout << "Got point " << result4.point.value() << std::endl;
        }
        else {
            std::cout << "Got no point." << std::endl;
        }
    }

    std::cout << "--- Intersection Tests Finished ---" << std::endl;
}
void run_segment_intersection_tests() {
    std::cout << "\n--- Running Segment Intersection Tests ---" << std::endl;
    using Status = geom::SegmentIntersectionResult2D<double>::Status;

    // Case 1: Simple intersection
    geom::Segment2d s1(geom::Point2d(0, 0), geom::Point2d(10, 10));
    geom::Segment2d s2(geom::Point2d(0, 10), geom::Point2d(10, 0));
    std::cout << "Test 1.1: Simple intersection... ";
    auto res1 = geom::intersection(s1, s2);
    if (res1.status == Status::INTERSECTING && res1.point.value() == geom::Point2d(5, 5)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Case 2: T-junction (endpoint on segment)
    geom::Segment2d s3(geom::Point2d(5, 5), geom::Point2d(15, 5));
    std::cout << "Test 1.2: T-junction intersection... ";
    auto res2 = geom::intersection(s1, s3);
    if (res2.status == Status::INTERSECTING && res2.point.value() == geom::Point2d(5, 5)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Case 3: Lines intersect, but segments do not
    geom::Segment2d s4(geom::Point2d(12, 12), geom::Point2d(20, 20));
    std::cout << "Test 1.3: No intersection (lines cross)... ";
    auto res3 = geom::intersection(s2, s4);
    if (res3.status == Status::NO_INTERSECTION) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Case 4: Collinear, but no overlap
    geom::Segment2d s5(geom::Point2d(0, 0), geom::Point2d(5, 0));
    geom::Segment2d s6(geom::Point2d(6, 0), geom::Point2d(10, 0));
    std::cout << "Test 1.4: Collinear with no overlap... ";
    auto res4 = geom::intersection(s5, s6);
    if (res4.status == Status::NO_INTERSECTION) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Case 5: Collinear with overlap
    geom::Segment2d s7(geom::Point2d(3, 0), geom::Point2d(8, 0));
    std::cout << "Test 1.5: Collinear with overlap... ";
    auto res5 = geom::intersection(s5, s7);
    geom::Segment2d expected_overlap(geom::Point2d(3, 0), geom::Point2d(5, 0));
    if (res5.status == Status::OVERLAPPING &&
        (res5.segment.value().p1() == expected_overlap.p1() && res5.segment.value().p2() == expected_overlap.p2())) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    std::cout << "--- Segment Intersection Tests Finished ---" << std::endl;
}

void run_polygon_tests() {
    std::cout << "\n--- Running Polygon Tests ---" << std::endl;

    // Case 1: Simple square
    std::vector<geom::Point2d> square_vertices = {
        {0, 0}, {5, 0}, {5, 5}, {0, 5}
    };
    geom::Polygon2d square(square_vertices);

    std::cout << "Test 1.1: Number of vertices... ";
    if (square.num_vertices() == 4) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    std::cout << "Test 1.2: Getting an edge... ";
    geom::Segment2d edge1 = square.edge(1); // The edge from (5,0) to (5,5)
    if (edge1.p1() == geom::Point2d(5, 0) && edge1.p2() == geom::Point2d(5, 5)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    std::cout << "Test 1.3: Getting the last edge (wraps around)... ";
    geom::Segment2d last_edge = square.edge(3); // The edge from (0,5) to (0,0)
    if (last_edge.p1() == geom::Point2d(0, 5) && last_edge.p2() == geom::Point2d(0, 0)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Test 2: Area calculation
    std::cout << "Test 2.1: Area of a 5x5 square... ";
    if (std::abs(square.area() - 25.0) < geom::Coord<double>::Epsilon) {
        std::cout << "SUCCESS. Area: " << square.area() << std::endl;
    }
    else { std::cout << "FAILED. Got " << square.area() << std::endl; }

    // Case 2: Triangle
    std::vector<geom::Point2d> triangle_vertices = { {0,0}, {10,0}, {5,5} };
    geom::Polygon2d triangle(triangle_vertices);
    // Area = 0.5 * base * height = 0.5 * 10 * 5 = 25
    std::cout << "Test 2.2: Area of a triangle... ";
    if (std::abs(triangle.area() - 25.0) < geom::Coord<double>::Epsilon) {
        std::cout << "SUCCESS. Area: " << triangle.area() << std::endl;
    }
    else { std::cout << "FAILED. Got " << triangle.area() << std::endl; }

    std::cout << "--- Polygon Tests Finished ---" << std::endl;
}

void run_convex_hull_tests() {
    std::cout << "\n--- Running Convex Hull Tests ---" << std::endl;

    std::vector<geom::Point2d> points = {
        {0, 0}, {5, 1}, {10, 0}, // Lower part
        {1, 5}, {9, 5},         // Upper part
        {5, 2}, {3, 3}, {7, 3}  // Internal points
    };

    auto hull_polygon_opt = geom::convex_hull(points);

    std::cout << "Test 1.1: Hull creation... ";
    if (hull_polygon_opt.has_value()) {
        std::cout << "SUCCESS" << std::endl;
    }
    else {
        std::cout << "FAILED to create hull." << std::endl;
        return;
    }

    auto hull_vertices = hull_polygon_opt.value().vertices();

    std::vector<geom::Point2d> expected_hull = {
        {0, 0}, {10, 0}, {9, 5}, {1, 5}
    };

    std::cout << "Test 1.2: Correct number of vertices... ";
    if (hull_vertices.size() == expected_hull.size()) {
        std::cout << "SUCCESS" << std::endl;
    }
    else {
        std::cout << "FAILED. Expected " << expected_hull.size() << ", got " << hull_vertices.size() << std::endl;
    }

    std::cout << "Test 1.3: Vertices match... ";
    bool match = true;
    for (const auto& expected_v : expected_hull) {
        bool found = false;
        for (const auto& actual_v : hull_vertices) {
            if (expected_v == actual_v) {
                found = true;
                break;
            }
        }
        if (!found) {
            match = false;
            break;
        }
    }

    if (match) {
        std::cout << "SUCCESS" << std::endl;
    }
    else {
        std::cout << "FAILED. Vertex sets do not match." << std::endl;
    }

    std::cout << "--- Convex Hull Tests Finished ---" << std::endl;
}
void run_point_in_polygon_tests() {
    std::cout << "\n--- Running Point in Polygon Tests ---" << std::endl;

    std::vector<geom::Point2d> vertices = { {0,0}, {10,0}, {10,10}, {5,15}, {0,10} };
    geom::Polygon2d poly(vertices);

    geom::Point2d p_inside(5, 5);
    geom::Point2d p_outside(15, 15);
    geom::Point2d p_on_edge(10, 5); // on the edge from (10,0) to (10,10)
    geom::Point2d p_on_vertex(0, 0);
    geom::Point2d p_tricky_case(6, 12); // near the concave vertex

    std::cout << "Test 1.1: Point inside... ";
    if (geom::contains(p_inside, poly)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 1.2: Point outside... ";
    if (!geom::contains(p_outside, poly)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 1.3: Point on edge... ";
    if (geom::contains(p_on_edge, poly)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 1.4: Point on vertex... ";
    if (geom::contains(p_on_vertex, poly)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 1.5: Tricky concave case (inside)... ";
    if (geom::contains(p_tricky_case, poly)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "--- Point in Polygon Tests Finished ---" << std::endl;
}


void run_ray_tests() {
    std::cout << "\n--- Running Ray Tests ---" << std::endl;

    geom::Ray2d ray = geom::Ray2d::from_points({ 1, 1 }, { 5, 4 });
    // Origin is (1,1). Direction is vector (4,3), normalized to (0.8, 0.6)

    // Test 1: Construction
    std::cout << "Test 1.1: Ray origin... ";
    if (ray.origin() == geom::Point2d(1, 1)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 1.2: Ray direction normalization... ";
    if (ray.direction() == geom::Vector2d(0.8, 0.6)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    // Test 2: Contains
    geom::Point2d p_on_ray = ray.point_at(10); // A point definitely on the ray
    geom::Point2d p_origin = ray.origin();
    geom::Point2d p_behind = ray.origin() - ray.direction(); // A point "behind" the origin
    geom::Point2d p_off_line(5, 5); // A point not on the supporting line

    std::cout << "Test 2.1: Contains point on ray... ";
    if (ray.contains(p_on_ray)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 2.2: Contains origin point... ";
    if (ray.contains(p_origin)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 2.3: Does NOT contain point behind origin... ";
    if (!ray.contains(p_behind)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "Test 2.4: Does NOT contain point off the line... ";
    if (!ray.contains(p_off_line)) std::cout << "SUCCESS" << std::endl;
    else std::cout << "FAILED" << std::endl;

    std::cout << "--- Ray Tests Finished ---" << std::endl;
}


void run_mixed_intersection_tests() {
    std::cout << "\n--- Running Mixed Intersection Tests ---" << std::endl;
    using Status = geom::SegmentIntersectionResult2D<double>::Status;

    geom::Line2d line = geom::Line2d::from_points({ 0, 5 }, { 10, 5 }); // Horizontal line y=5

    // Case 1: Line intersects segment
    geom::Segment2d seg1({ 5, 0 }, { 5, 10 });
    std::cout << "Test 1.1: Line intersects segment... ";
    auto res1 = geom::intersection(line, seg1);
    if (res1.status == Status::INTERSECTING && res1.point.value() == geom::Point2d(5, 5)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Case 2: Line does not intersect segment (but supporting line does)
    geom::Segment2d seg2({ 5, 6 }, { 5, 10 });
    std::cout << "Test 1.2: Line misses segment... ";
    auto res2 = geom::intersection(line, seg2);
    if (res2.status == Status::NO_INTERSECTION) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // Case 3: Segment is coincident with line
    geom::Segment2d seg3({ -2, 5 }, { 4, 5 });
    std::cout << "Test 1.3: Segment coincident with line... ";
    auto res3 = geom::intersection(line, seg3);
    if (res3.status == Status::OVERLAPPING &&
        res3.segment.value().p1() == seg3.p1() && res3.segment.value().p2() == seg3.p2()) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }


    std::cout << "--- Mixed Intersection Tests Finished ---" << std::endl;

    // ... run_mixed_intersection_tests() ֆունկցիայի ներսում ...
// ... seg3-ի թեստից հետո ...

     // --- Line-Ray Tests ---
    // *** ՈՒՂՂՎԱԾ Է ***
    geom::Ray2d ray = geom::Ray2d::from_point_direction({ 5, 0 }, { 0, 1 }); // Vertical ray starting at (5,0) going upwards

    std::cout << "Test 2.1: Line intersects ray... ";
    auto res4 = geom::intersection(line, ray); // line is y=5
    if (res4.status == Status::INTERSECTING && res4.point.value() == geom::Point2d(5, 5)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // *** ՈՒՂՂՎԱԾ Է ***
    std::cout << "Test 2.2: Line misses ray (intersects behind)... ";
    geom::Ray2d ray2 = geom::Ray2d::from_point_direction({ 5, 6 }, { 0, 1 }); // Starts above the line
    auto res5 = geom::intersection(line, ray2);
    if (res5.status == Status::NO_INTERSECTION) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    // *** ՈՒՂՂՎԱԾ Է ***
    std::cout << "Test 2.3: Ray coincident with line... ";
    geom::Ray2d ray3 = geom::Ray2d::from_point_direction({ -2, 5 }, { 1, 0 }); // Coincident, goes right
    auto res6 = geom::intersection(line, ray3);
    if (res6.status == Status::OVERLAPPING) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }


    std::cout << "--- Mixed Intersection Tests Finished ---" << std::endl;

    // ... run_mixed_intersection_tests() ֆունկցիայի ներսում ...
// ... ray3-ի թեստից հետո ...

    std::cout << std::endl; // Add a separator

    // --- Segment-Ray Tests ---
    geom::Segment2d seg({ -5, 0 }, { 5, 10 });

    std::cout << "Test 3.1: Segment intersects ray... ";
    geom::Ray2d ray_s1 = geom::Ray2d::from_point_direction({ -10, 5 }, { 1, 0 }); // Horizontal ray from left
    auto res7 = geom::intersection(seg, ray_s1);
    if (res7.status == Status::INTERSECTING && res7.point.value() == geom::Point2d(0, 5)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    std::cout << "Test 3.2: Segment misses ray (ray starts after)... ";
    geom::Ray2d ray_s2 = geom::Ray2d::from_point_direction({ 1, 5 }, { 1, 0 }); // Horizontal ray starts after intersection
    auto res8 = geom::intersection(seg, ray_s2);
    if (res8.status == Status::NO_INTERSECTION) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }

    std::cout << "Test 3.3: Collinear overlap (ray starts inside segment)... ";
    geom::Segment2d seg_col({ 0, 0 }, { 10, 0 });
    geom::Ray2d ray_col = geom::Ray2d::from_point_direction({ 5, 0 }, { 1, 0 });
    auto res9 = geom::intersection(seg_col, ray_col);
    if (res9.status == Status::OVERLAPPING && res9.segment.value().p1() == geom::Point2d(5, 0) && res9.segment.value().p2() == geom::Point2d(10, 0)) {
        std::cout << "SUCCESS" << std::endl;
    }
    else { std::cout << "FAILED" << std::endl; }


    // ... ֆունկցիայի վերջը ...
}



// ... main() ֆունկցիայի մեջ ...
int main() {
    run_vector_tests();
    run_line_tests();
    run_segment_tests();
    run_intersection_tests(); // Ավելացնում ենք կանչը
    run_segment_intersection_tests();
    run_polygon_tests();
    run_convex_hull_tests();
    run_point_in_polygon_tests();
    run_ray_tests();
    run_mixed_intersection_tests();
    return 0;
}