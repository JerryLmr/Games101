#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>
#include <iostream>
#include <string>

struct RunStats
{
    double setupMs = 0.0;
    double renderMs = 0.0;
    double totalMs = 0.0;
};

static const char* splitMethodName(BVHAccel::SplitMethod splitMethod)
{
    return splitMethod == BVHAccel::SplitMethod::SAH ? "SAH" : "NAIVE";
}

static void printUsage(const char* programName)
{
    std::cout << "Usage:\n"
              << "  " << programName << "\n"
              << "  " << programName << " --bvh sah\n"
              << "  " << programName << " --bvh naive\n"
              << "  " << programName << " --compare\n";
}

static bool parseSplitMethod(const std::string& value,
                             BVHAccel::SplitMethod& splitMethod)
{
    if (value == "sah") {
        splitMethod = BVHAccel::SplitMethod::SAH;
        return true;
    }
    if (value == "naive") {
        splitMethod = BVHAccel::SplitMethod::NAIVE;
        return true;
    }
    return false;
}

static RunStats runOnce(BVHAccel::SplitMethod splitMethod,
                        const std::string& label,
                        const std::string& outputFile)
{
    std::cout << "\n=== " << label << " (" << splitMethodName(splitMethod)
              << ") ===\n";

    auto totalStart = std::chrono::steady_clock::now();
    auto setupStart = totalStart;

    Scene scene(1280, 960);

    MeshTriangle bunny("../models/bunny/bunny.obj", splitMethod);

    scene.Add(&bunny);
    scene.Add(std::make_unique<Light>(Vector3f(-20, 70, 20), 1));
    scene.Add(std::make_unique<Light>(Vector3f(20, 70, 20), 1));
    scene.buildBVH(splitMethod);

    auto setupStop = std::chrono::steady_clock::now();

    Renderer r;

    auto renderStart = std::chrono::steady_clock::now();
    r.Render(scene, outputFile);
    auto renderStop = std::chrono::steady_clock::now();

    auto totalStop = renderStop;

    RunStats stats;
    stats.setupMs =
        std::chrono::duration<double, std::milli>(setupStop - setupStart).count();
    stats.renderMs =
        std::chrono::duration<double, std::milli>(renderStop - renderStart).count();
    stats.totalMs =
        std::chrono::duration<double, std::milli>(totalStop - totalStart).count();

    std::cout << "\nRender complete (" << label << "):\n";
    std::cout << "  Output file      : " << outputFile << "\n";
    std::cout << "  Scene setup time : " << stats.setupMs << " ms\n";
    std::cout << "  Render time      : " << stats.renderMs << " ms\n";
    std::cout << "  Total time       : " << stats.totalMs << " ms\n";

    return stats;
}

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{
    BVHAccel::SplitMethod splitMethod = BVHAccel::SplitMethod::SAH;

    if (argc == 1) {
        runOnce(splitMethod, "default", "binary.ppm");
        return 0;
    }

    if (argc == 2 && std::string(argv[1]) == "--compare") {
        RunStats naiveStats =
            runOnce(BVHAccel::SplitMethod::NAIVE, "naive", "binary_naive.ppm");
        RunStats sahStats =
            runOnce(BVHAccel::SplitMethod::SAH, "sah", "binary_sah.ppm");

        std::cout << "\n=== Comparison ===\n";
        std::cout << "NAIVE setup : " << naiveStats.setupMs << " ms\n";
        std::cout << "SAH setup   : " << sahStats.setupMs << " ms\n";
        std::cout << "NAIVE render: " << naiveStats.renderMs << " ms\n";
        std::cout << "SAH render  : " << sahStats.renderMs << " ms\n";
        std::cout << "NAIVE total : " << naiveStats.totalMs << " ms\n";
        std::cout << "SAH total   : " << sahStats.totalMs << " ms\n";
        return 0;
    }

    if (argc == 3 && std::string(argv[1]) == "--bvh" &&
        parseSplitMethod(argv[2], splitMethod)) {
        runOnce(splitMethod, argv[2], "binary.ppm");
        return 0;
    }

    printUsage(argv[0]);

    return 1;
}
