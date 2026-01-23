#include "NavigationSystem.h"
#include "Triangulation.h"
#include "AStarPathfinder.h"
#include "FunnelAlgorithm.h"
#include <iostream>

using namespace Geom;

NavigationSystem::NavigationSystem()
    : window(sf::VideoMode(1200, 800), "Autonomous Navigation System"),
    mode(Mode::DRAWING_OBSTACLES),
    drawingObstacle(false),
    draggingObstacleIndex(-1),
    navMeshBuilt(false),
    hasStart(false),
    hasGoal(false),
    isPanning(false),
    zoomLevel(1.0f),
    showNavMesh(true),
    showPath(true),
    showBounds(false),
    obstacleColor(sf::Color::Red),
    navMeshColor(sf::Color(100, 100, 255, 150)),
    pathColor(sf::Color(0, 255, 0)),
    startColor(sf::Color::Blue),
    goalColor(sf::Color::Magenta),
    boundsColor(sf::Color(150, 150, 150))
{
    bounds.addVertex(Point<double>(50, 50));
    bounds.addVertex(Point<double>(1150, 50));
    bounds.addVertex(Point<double>(1150, 750));
    bounds.addVertex(Point<double>(50, 750));

    view.setSize(1200, 800);
    view.setCenter(600, 400);
    window.setView(view);
    window.setFramerateLimit(60);
}

void NavigationSystem::run() {
    while (window.isOpen()) {
        sf::Time dt = deltaClock.restart();
        double deltaTime = dt.asSeconds();
        if (deltaTime > 0.1) deltaTime = 0.1;

        handleEvents();
        agent.update(deltaTime);
        render();
    }
}

void NavigationSystem::buildNavMesh() {
    std::vector<Triangle<double>> walkableTriangles =
        triangulateWithHoles(bounds, obstacles);

    navMesh.build(walkableTriangles);
    navMeshBuilt = true;

    if (hasGoal) {
        findPath();
    }
}

void NavigationSystem::findPath() {
    if (!navMeshBuilt || !hasGoal) return;

    Point<double> searchStart = startPoint;
    if (agent.isActive()) {
        searchStart = agent.getPosition();
    }
    else if (!hasStart) {
        return;
    }

    bool startIsValid = true;
    for (const auto& obs : obstacles) {
        if (isPointInPolygon(obs, searchStart)) {
            startIsValid = false;
            break;
        }
    }

    if (!startIsValid) {
        agent.stop();
        smoothPath.clear();
        return;
    }

    AStarPathfinder<double> pathfinder(&navMesh);
    trianglePath = pathfinder.findPath(searchStart, goalPoint);

    if (trianglePath.empty()) {
        smoothPath.clear();
        agent.stop();
        return;
    }

    FunnelAlgorithm<double> funnel(&navMesh);
    smoothPath = funnel.smoothPath(trianglePath, searchStart, goalPoint);

    agent.setPath(smoothPath);
}

void NavigationSystem::handleEvents() {
    sf::Event event;
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed) window.close();

        if (event.type == sf::Event::MouseWheelScrolled) {
            if (event.mouseWheelScroll.delta > 0) { view.zoom(0.9f); zoomLevel *= 0.9f; }
            else { view.zoom(1.1f); zoomLevel *= 1.1f; }
            window.setView(view);
        }

        if (event.type == sf::Event::KeyPressed) {
            if (event.key.code == sf::Keyboard::Escape) window.close();

            if (event.key.code == sf::Keyboard::Left) { view.move(-50.f * zoomLevel, 0); window.setView(view); }
            if (event.key.code == sf::Keyboard::Right) { view.move(50.f * zoomLevel, 0); window.setView(view); }
            if (event.key.code == sf::Keyboard::Up) { view.move(0, -50.f * zoomLevel); window.setView(view); }
            if (event.key.code == sf::Keyboard::Down) { view.move(0, 50.f * zoomLevel); window.setView(view); }

            if (event.key.code == sf::Keyboard::D) { mode = Mode::DRAWING_OBSTACLES; std::cout << "Mode: Drawing" << std::endl; }
            if (event.key.code == sf::Keyboard::B) buildNavMesh();
            if (event.key.code == sf::Keyboard::S) { mode = Mode::SELECTING_START; std::cout << "Mode: Start" << std::endl; }
            if (event.key.code == sf::Keyboard::G) { mode = Mode::SELECTING_GOAL; std::cout << "Mode: Goal" << std::endl; }

            if (event.key.code == sf::Keyboard::C) {
                obstacles.clear();
                currentObstacle = Polygon<double>();
                navMeshBuilt = false; hasStart = false; hasGoal = false;
                smoothPath.clear();
                agent.stop();
                std::cout << "Cleared" << std::endl;
            }

            if (event.key.code == sf::Keyboard::R) {
                view.setSize(1200, 800); view.setCenter(600, 400); zoomLevel = 1.0f; window.setView(view);
            }
            if (event.key.code == sf::Keyboard::M) showNavMesh = !showNavMesh;
            if (event.key.code == sf::Keyboard::P) showPath = !showPath;
            if (event.key.code == sf::Keyboard::O) showBounds = !showBounds;
        }

        if (event.type == sf::Event::MouseButtonPressed) {
            if (event.mouseButton.button == sf::Mouse::Middle) {
                isPanning = true;
                lastMousePos = sf::Vector2i(event.mouseButton.x, event.mouseButton.y);
            }

            if (event.mouseButton.button == sf::Mouse::Left) {
                sf::Vector2f worldPos = window.mapPixelToCoords(sf::Vector2i(event.mouseButton.x, event.mouseButton.y));
                Point<double> clickPos(worldPos.x, worldPos.y);
                lastMousePos = sf::Vector2i(event.mouseButton.x, event.mouseButton.y);

                bool clickedObstacle = false;
                for (int i = obstacles.size() - 1; i >= 0; i--) {
                    if (isPointInPolygon(obstacles[i], clickPos)) {
                        draggingObstacleIndex = i;
                        mode = Mode::DRAGGING_OBSTACLE;
                        clickedObstacle = true;
                        break;
                    }
                }

                if (!clickedObstacle) {
                    if (mode == Mode::DRAWING_OBSTACLES) {
                        currentObstacle.addVertex(clickPos);
                        drawingObstacle = true;
                    }
                    else if (mode == Mode::SELECTING_START) {
                        startPoint = clickPos;
                        hasStart = true;
                        agent.setPosition(startPoint);
                        agent.stop();
                        if (hasGoal && navMeshBuilt) findPath();
                    }
                    else if (mode == Mode::SELECTING_GOAL) {
                        goalPoint = clickPos;
                        hasGoal = true;
                        if (navMeshBuilt) findPath();
                    }
                }
            }

            if (event.mouseButton.button == sf::Mouse::Right) {
                if (drawingObstacle) {
                    if (currentObstacle.size() >= 3) {
                        obstacles.push_back(currentObstacle);
                        std::cout << "Obstacle added" << std::endl;
                    }
                    currentObstacle = Polygon<double>();
                    drawingObstacle = false;
                }
            }
        }

        if (event.type == sf::Event::MouseButtonReleased) {
            if (event.mouseButton.button == sf::Mouse::Middle) isPanning = false;

            if (event.mouseButton.button == sf::Mouse::Left) {
                if (draggingObstacleIndex != -1) {
                    draggingObstacleIndex = -1;
                    buildNavMesh();
                }
            }
        }

        if (event.type == sf::Event::MouseMoved) {
            sf::Vector2i currentMousePos(event.mouseMove.x, event.mouseMove.y);
            sf::Vector2f currWorld = window.mapPixelToCoords(currentMousePos);
            sf::Vector2f lastWorld = window.mapPixelToCoords(lastMousePos);
            sf::Vector2f delta = currWorld - lastWorld;

            if (draggingObstacleIndex != -1) {
                for (auto& vertex : obstacles[draggingObstacleIndex].vertices) {
                    vertex.x += delta.x;
                    vertex.y += delta.y;
                }
            }
            else if (isPanning) {
                sf::Vector2i pixelDelta = lastMousePos - currentMousePos;
                view.move(pixelDelta.x * zoomLevel, pixelDelta.y * zoomLevel);
                window.setView(view);
            }

            lastMousePos = currentMousePos;
        }
    }
}

void NavigationSystem::render() {
    window.clear(sf::Color::White);

    if (showBounds) drawPolygon(bounds, boundsColor);
    for (const auto& obs : obstacles) drawPolygon(obs, obstacleColor, true);

    if (drawingObstacle && currentObstacle.size() > 0) {
        drawPolygon(currentObstacle, sf::Color(255, 150, 150));
        for (auto& p : currentObstacle.vertices) drawPoint(p, sf::Color::Red, 5.f);
    }

    if (navMeshBuilt && showNavMesh) {
        for (const auto& navTri : navMesh.getTriangles()) drawTriangle(navTri.tri, navMeshColor);
    }

    if (hasStart) drawPoint(startPoint, startColor, 8.f);
    if (hasGoal) drawPoint(goalPoint, goalColor, 8.f);

    if (!smoothPath.empty() && showPath) {
        for (size_t i = 0; i < smoothPath.size() - 1; i++) {
            drawLine(smoothPath[i], smoothPath[i + 1], pathColor, 4.f);
        }
        for (const auto& p : smoothPath) drawPoint(p, pathColor, 5.f);
    }

    if (agent.isActive() || hasStart) {
        if (agent.isActive()) drawPoint(agent.getPosition(), sf::Color::Black, 10.f);
    }

    window.display();
}

// Helpers
void NavigationSystem::drawPoint(const Point<double>& p, const sf::Color& color, float radius) {
    sf::CircleShape circle(radius);
    circle.setFillColor(color);
    circle.setPosition(static_cast<float>(p.x - radius), static_cast<float>(p.y - radius));
    window.draw(circle);
}

void NavigationSystem::drawLine(const Point<double>& p1, const Point<double>& p2, const sf::Color& color, float thickness) {
    sf::Vertex line[] = {
        sf::Vertex(sf::Vector2f(static_cast<float>(p1.x), static_cast<float>(p1.y)), color),
        sf::Vertex(sf::Vector2f(static_cast<float>(p2.x), static_cast<float>(p2.y)), color)
    };
    window.draw(line, 2, sf::Lines);
}

void NavigationSystem::drawPolygon(const Polygon<double>& poly, const sf::Color& color, bool filled) {
    if (poly.size() < 2) return;
    if (filled) {
        sf::ConvexShape shape;
        shape.setPointCount(poly.size());
        for (size_t i = 0; i < poly.size(); i++) {
            shape.setPoint(i, sf::Vector2f(static_cast<float>(poly.vertices[i].x), static_cast<float>(poly.vertices[i].y)));
        }
        shape.setFillColor(sf::Color(color.r, color.g, color.b, 100));
        shape.setOutlineColor(color);
        shape.setOutlineThickness(2.f);
        window.draw(shape);
    }
    else {
        for (size_t i = 0; i < poly.size(); i++) {
            size_t next = (i + 1) % poly.size();
            drawLine(poly.vertices[i], poly.vertices[next], color);
        }
    }
}

void NavigationSystem::drawTriangle(const Triangle<double>& tri, const sf::Color& color) {
    drawLine(tri.a, tri.b, color, 1.f);
    drawLine(tri.b, tri.c, color, 1.f);
    drawLine(tri.c, tri.a, color, 1.f);
}