#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>
#include <variant>
#include <functional>

// Forward declaration.
struct Scene;
struct TrianglePrim {
    size_t meshIdx;
    size_t v1, v2, v3;
};

struct SpherePrim {
    size_t sphereIdx;
};

struct Primitive {
    std::variant<TrianglePrim, SpherePrim> p;
    glm::vec3 center;
};

struct SAHCuts {
    std::array<std::vector<AxisAlignedBox>, 3> cuts;
    size_t chosenDim;
    size_t chosenInd;
};

struct Pair{
    int node;
    float t;
};

struct Node {
    AxisAlignedBox aabb;

    // First element: 0 -> internal node, 1 -> leaf node (index 0)
    // If internal node: depth, beg, end, left, right (index 1-5)
    // If leaf node: depth, beg, end (index 1-3)

    // beg, end - index-pointers into the primitives vector (beg - incl, end - excl)
    // left, right - indexes of child nodes
    std::vector<size_t> data;
};

const size_t MAX_DEPTH = 16;
const size_t NUM_OF_BINS = 8;

const std::array<std::function<bool(const Primitive&, const Primitive&)>, 3> comparators = {
    [](const auto& a, const auto& b) { return a.center.x < b.center.x; },
    [](const auto& a, const auto& b) { return a.center.y < b.center.y; },
    [](const auto& a, const auto& b) { return a.center.z < b.center.z; },
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Set recursion level
    void setRecursionLevel(int level, bool debug);

    void debugDrawSAHSplits(int level, int axis);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    int m_recursionLevel;
    int RECURSION_LEVEL;


private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Primitive> primitives;
    std::vector<Node> nodes;
    size_t root;

    std::vector<std::vector<SAHCuts>> sahCutsPerLevel;

    // this function should split the triangles in the given range, returns the split index
    std::function<
        size_t(std::vector<Primitive>& prims, size_t beg, size_t end, size_t depth, Scene* scene, std::vector<SAHCuts>& debugCuts)
        > splitFunc;

    size_t createBVH(size_t beg, size_t end, size_t depth);

};