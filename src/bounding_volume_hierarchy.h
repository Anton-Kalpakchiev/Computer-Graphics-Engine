#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>
#include <variant>

// Forward declaration.
struct Scene;
typedef glm::uvec3 Triangle;

struct PrimitiveVariant {
    std::variant<Triangle, Sphere> primitive;
    glm::vec3 center;
};

struct Internal {
    AxisAlignedBox aabb;
    size_t left, right;    
};

struct Node {
    std::variant<Internal, PrimitiveVariant> val;
    size_t level; // used for debug
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Vertex> vertices;
    std::vector<PrimitiveVariant> primitives;
    std::vector<Node> nodes;
    size_t root;

    size_t createBVH(size_t beg, size_t end, size_t splitBy, size_t depth);
};