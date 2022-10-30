#include "bvh_interface.h"
#include "bounding_volume_hierarchy.h"

//! DON'T TOUCH THIS FILE!

BvhInterface::BvhInterface(Scene* pScene, const Features& features)
{
    m_impl = new BoundingVolumeHierarchy(pScene, features);
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BvhInterface::numLevels() const
{
    return m_impl->numLevels();
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BvhInterface::numLeaves() const
{
    return m_impl->numLeaves();
}


// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BvhInterface::debugDrawLevel(int level)
{
    m_impl->debugDrawLevel(level);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function 
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BvhInterface::debugDrawLeaf(int leafIdx)
{
    m_impl->debugDrawLeaf(leafIdx);
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BvhInterface::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    return m_impl->intersect(ray, hitInfo, features);
}
