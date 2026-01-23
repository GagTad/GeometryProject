#pragma once
#include <SFML/Graphics.hpp>
#include "geometry.h"
#include "NavMesh.h"
#include "Agent.h"
#include <vector>

class NavigationSystem {
public:
    NavigationSystem();
    void run();

private:

    enum class Mode {
        DRAWING_OBSTACLES,
        SELECTING_START,
        SELECTING_GOAL,
        SHOWING_PATH,
        DRAGGING_OBSTACLE
    };

    void handleEvents();
    void render();
    void buildNavMesh();
    void findPath();

    void drawPoint(const Geom::Point<double>& p, const sf::Color& color, float radius = 5.f);
    void drawLine(const Geom::Point<double>& p1, const Geom::Point<double>& p2, const sf::Color& color, float thickness = 2.f);
    void drawPolygon(const Geom::Polygon<double>& poly, const sf::Color& color, bool filled = false);
    void drawTriangle(const Geom::Triangle<double>& tri, const sf::Color& color);

    sf::RenderWindow window;
    sf::View view;
    Mode mode;
    sf::Clock deltaClock;

    Geom::Polygon<double> bounds;
    std::vector<Geom::Polygon<double>> obstacles;
    Geom::Polygon<double> currentObstacle;
    bool drawingObstacle;
    int draggingObstacleIndex;

    Geom::NavMesh<double> navMesh;
    bool navMeshBuilt;

    Geom::Point<double> startPoint;
    Geom::Point<double> goalPoint;
    bool hasStart;
    bool hasGoal;
    std::vector<Geom::Point<double>> smoothPath;
    std::vector<int> trianglePath;

    Geom::Agent agent;

    bool isPanning;
    sf::Vector2i lastMousePos;
    float zoomLevel;

    bool showNavMesh;
    bool showPath;
    bool showBounds;

    sf::Color obstacleColor;
    sf::Color navMeshColor;
    sf::Color pathColor;
    sf::Color startColor;
    sf::Color goalColor;
    sf::Color boundsColor;
};