#include "NavigationSystem.h"
#include <iostream>

int main() {
    std::cout << "=== AUTONOMOUS NAVIGATION SYSTEM (DYNAMIC) ===" << std::endl;
    std::cout << "Workflow:" << std::endl;
    std::cout << "  1. Draw Obstacles [D] or DRAG existing ones with mouse!" << std::endl;
    std::cout << "  2. Press [B] to build NavMesh initially" << std::endl;
    std::cout << "  3. Set [S]tart and [G]oal -> Agent moves" << std::endl;
    std::cout << "  4. Move obstacles while agent moves -> Path updates automatically!" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "  D/B/S/G - Drawing/Build/Start/Goal" << std::endl;
    std::cout << "  C - Clear all" << std::endl;
    std::cout << "  R - Reset view" << std::endl;
    std::cout << "==============================================" << std::endl;

    NavigationSystem app;
    app.run();

    return 0;
}
