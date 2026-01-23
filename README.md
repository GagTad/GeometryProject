#  Autonomous Navigation System

A high-performance **2D Pathfinding Engine** built from scratch using **C++17** and **SFML**. This project implements advanced computational geometry algorithms to generate navigation meshes, calculate optimal paths, and handle dynamic environments in real-time.

![C++](https://img.shields.io/badge/Language-C%2B%2B17-blue)
![SFML](https://img.shields.io/badge/Library-SFML-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![License](https://img.shields.io/badge/License-MIT-orange)


<img width="762" height="411" alt="image" src="https://github.com/user-attachments/assets/26867909-3c0b-46ec-be23-625e7130546f" />

##  Overview

Unlike simple tile-based pathfinders, this system uses a **Navigation Mesh (NavMesh)** approach. It converts the walkable environment into a graph of convex polygons (triangles), allowing for smoother, more natural movement and significantly lower memory usage in large open spaces.

The engine features a **Dynamic World** capability: obstacles can be interactively moved, triggering immediate mesh regeneration and path recalculation without stalling the simulation.

##  Key Features

*   **Custom Geometry Engine:** Implemented a robust math library from scratch (Vectors, Points, Polygons) without relying on heavy external physics engines.
*   **NavMesh Generation:**
    *   **Grid-Based Triangulation:** Discretizes space and triangulates walkable areas.
    *   **Hole Detection:** Automatically identifies and removes triangles intersecting with obstacles.
    *   **AABB Optimization:** Uses Axis-Aligned Bounding Boxes for fast rejection during collision checks.
*   **Pathfinding Algorithms:**
    *   **A* (A-Star):** Finds the optimal sequence of triangles (corridor) to reach the goal.
    *   **Funnel Algorithm (SSFA):** Post-processes the A* output using a "string pulling" technique to find the shortest physical path within the corridor.
*   **Dynamic Environment:**
    *   **Interactive Editing:** Drag & Drop obstacles in real-time.
    *   **Live Re-pathing:** Agents intelligently recalculate paths from their current velocity and position when the world changes.

##  Tech Stack & Architecture

The project follows a modular **Clean Architecture** to ensure maintainability:

*   **Core:** `NavigationSystem` (Manages the game loop, input, and rendering).
*   **Geometry:** `Point.h`, `Vector.h`, `Polygon.h` (Template-based math primitives).
*   **Algorithms:** `Triangulation.h`, `AStarPathfinder.h`, `FunnelAlgorithm.h` (Pure logic, separated from rendering).
*   **Entities:** `Agent.h` (Handles physics-based movement and state).

##  Controls

| Key / Action | Description |
| :--- | :--- |
| **`D`** | Enter **Draw Mode** (Left click to add points, Right click to finish obstacle) |
| **`S`** | Set **Start** Point |
| **`G`** | Set **Goal** Point |
| **`B`** | **Build** / Rebake NavMesh manually |
| **`C`** | **Clear** All |
| **`M`** | Toggle **NavMesh** Visualization |
| **`P`** | Toggle **Path** Visualization |
| **Left Click + Drag** | **Move Obstacles** (Dynamic updates) |
| **Mouse Wheel** | Zoom In/Out |
| **Arrows** | Pan Camera |

##  Getting Started

### Prerequisites
*   C++ Compiler (supporting C++17)
*   SFML 2.5+ Library

### Installation (Visual Studio)
1.  Clone the repository.
2.  Open the solution in Visual Studio.
3.  Ensure SFML include and library paths are configured.
4.  Build and Run.

### Installation (CLI / Linux)
```bash
# Install SFML
sudo apt-get install libsfml-dev

# Compile
g++ -std=c++17 src/*.cpp -o nav_system -lsfml-graphics -lsfml-window -lsfml-system

# Run
./nav_system
