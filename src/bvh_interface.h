#pragma once
#include "config.h"
#include <array>
#include <span>

//! DON'T TOUCH THIS FILE! !//

// Forward declaration.
class BoundingVolumeHierarchy;
struct Scene;

class BvhInterface {
public:

    // Constructor. Receives the scene and builds the bounding volume hierarchy
    BvhInterface(Scene* pScene, const Features& features);


    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;


    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;


    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Set recursion level
    void setRecursionLevel(int level) const;

    // Set debug recursion level
    void setDebugRecursionLevel(int level) const;

    // Visual Debug SAH + binning: draw considered splits and the chosen one per level
    void debugDrawSAHSplits(int level, int axis);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

private:
    BoundingVolumeHierarchy* m_impl;
};
