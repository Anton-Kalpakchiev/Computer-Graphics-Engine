#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <algorithm>
#include <optional>
#include <span>
#include <functional>
#include <ranges>
#include <framework/variant_helper.h>
#include <variant>
#include <optional>
#include <stack>
#include <iostream>

glm::vec3 triangleCenter(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
    return a + b + c / 3.f;
}

std::optional<AxisAlignedBox> getBoundingBox(std::span<Primitive> primitives, size_t beg, size_t end) {
    std::optional<AxisAlignedBox> res;
    
    for (auto it = primitives.begin() + beg; it != primitives.begin() + end; it++) {
        float xMin, xMax, yMin, yMax, zMin, zMax;

        std::visit(make_visitor(
            [&](const TrianglePrim& t) {
                const auto& p1 = t.v1->position;
                const auto& p2 = t.v2->position;
                const auto& p3 = t.v3->position;
                xMin = std::min({ p1.x, p2.x, p3.x });
                xMax = std::max({ p1.x, p2.x, p3.x });
                yMin = std::min({ p1.y, p2.y, p3.y });
                yMax = std::max({ p1.y, p2.y, p3.y });
                zMin = std::min({ p1.z, p2.z, p3.z });
                zMax = std::max({ p1.z, p2.z, p3.z });
            },
            [&](const SpherePrim& s) {
                xMin = s.c->x - s.r;
                xMax = s.c->x + s.r;
                yMin = s.c->y - s.r;
                yMax = s.c->y + s.r;
                zMin = s.c->z - s.r;
                zMax = s.c->z + s.r;
            }
        ), it->p);

        auto v = res.value_or(AxisAlignedBox { glm::vec3(xMin, yMin, zMin), glm::vec3(xMax, yMax, zMax) });
        v.lower.x = std::min(v.lower.x, xMin);
        v.lower.y = std::min(v.lower.y, yMin);
        v.lower.z = std::min(v.lower.z, zMin);
        v.upper.x = std::max(v.upper.x, xMax);
        v.upper.y = std::max(v.upper.y, yMax);
        v.upper.z = std::max(v.upper.z, zMax);
        res = v;
    }
    return res;
}

size_t BoundingVolumeHierarchy::createBVH(size_t beg, size_t end, size_t splitBy, size_t depth) {
    m_numLevels = std::max(m_numLevels, (int)depth + 1);
    auto aabb = getBoundingBox(primitives, beg, end).value();
    if (depth + 1 == MAX_DEPTH || beg + 1 == end) {
        nodes.push_back(Node { aabb, { true, depth, beg, end } });
        m_numLeaves++;
        return nodes.size() - 1;
    }
    auto byX = [](const auto& a, const auto& b) { return a.center.x < b.center.x; };
    auto byY = [](const auto& a, const auto& b) { return a.center.y < b.center.y; };
    auto byZ = [](const auto& a, const auto& b) { return a.center.z < b.center.z; };
    const std::function<bool(const Primitive&, const Primitive&)> comparators[] = {byX, byY, byZ};

    size_t mid = beg + (end - beg) / 2;
    std::nth_element(primitives.begin() + beg, primitives.begin() + mid, primitives.begin() + end, comparators[splitBy]);

    auto left = createBVH(beg, mid, (splitBy + 1) % 3, depth + 1);
    auto right = createBVH(mid, end, (splitBy + 1) % 3, depth + 1);
    nodes.push_back(Node {aabb, { false, depth, beg, end, left, right } });
    return nodes.size() - 1;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    
    m_numLevels = 0;

    // Get all triangles of the scene into the BVH
    for (auto& mesh : pScene->meshes) {

        std::transform(mesh.triangles.begin(), mesh.triangles.end(), std::back_inserter(primitives), [&](const auto& t) {
            auto p1 = &mesh.vertices[t.x];
            auto p2 = &mesh.vertices[t.y];
            auto p3 = &mesh.vertices[t.z];
            return Primitive { TrianglePrim { p1, p2, p3 }, &mesh.material, triangleCenter(p1->position, p2->position, p3->position) };
        });
    }
    // Get all spheres of the scene into the BVH
    std::transform(pScene->spheres.begin(), pScene->spheres.end(), std::back_inserter(primitives), [](auto& s) {
        return Primitive { SpherePrim { &s.center, s.radius }, &s.material, glm::vec3(s.center) };
    });

    m_numLeaves = 0;
    m_numLevels = 0;

    // We have all the primitives and their centers in the primitves vector
    // Create the BVH itself
    root = createBVH(0, primitives.size(), 0, 0);
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    for (const auto& node : nodes) {
        if (node.data[1] == level) {
            drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(0.9f));
        }
    }
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
    for (const auto& node : nodes) {
        if (node.data[0]) leafIdx--;
        if (leafIdx == 0) {
            drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(.2f, 1.f, .2f), .5f);
            setColor({1.f, 1.f, 0.f});
            auto beg = primitives.begin() + node.data[2];
            auto end = primitives.begin() + node.data[3];
            for (const auto& p : std::ranges::subrange(beg, end)) {
                std::visit(make_visitor(
                    [&](const TrianglePrim& t) {
                        drawTriangle(*t.v1, *t.v2, *t.v3);
                    },
                    [&](const SpherePrim& s) {
                        drawSphere(Sphere { *s.c, s.r, *p.mat });
                    }
                ), p.p);
            }
            break;
        }
    }
}

std::optional<Primitive> getIntersecting(auto beg, auto end, Ray& ray, HitInfo& hitInfo) {
    std::optional<Primitive> res;
    auto range = std::ranges::subrange(beg, end);
    for (const auto& prim : range) {
        auto hit =  std::visit(make_visitor(
            [&](const TrianglePrim& t) {
                return intersectRayWithTriangle(t.v1->position, t.v2->position, t.v3->position, ray, hitInfo);
            },
            [&](const SpherePrim& s) {
                return intersectRayWithShape(Sphere { *s.c, s.r }, ray, hitInfo);
            }
        ), prim.p);

        if (hit) {
            res = prim;
        }
    }
    return res;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    std::optional<Primitive> prim;

    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        
        prim = getIntersecting(primitives.begin(), primitives.end(), ray, hitInfo);
        
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.

        auto stack = std::stack<Node>();
        stack.push(nodes[root]);

        while(!stack.empty()){

            auto parent = stack.top();
            stack.pop();

            if (parent.data[0] == 1) {

                size_t beg = parent.data[2];
                size_t end = parent.data[3];
                auto maybePrim = getIntersecting(primitives.begin() + beg, primitives.begin() + end, ray, hitInfo);
                if (maybePrim.has_value()) {
                    prim = maybePrim.value();
                }

            } else {

                auto left = nodes[parent.data[4]];
                auto right = nodes[parent.data[5]];

                auto rollBack = ray.t;

                ray.t = std::numeric_limits<float>::max();
                bool leftBox = intersectRayWithShape(left.aabb, ray);

                ray.t = std::numeric_limits<float>::max();
                bool rightBox = intersectRayWithShape(right.aabb, ray);

                ray.t = rollBack;
                
                if (leftBox) stack.push(left);
                if (rightBox) stack.push(right);
            }
        }
    }
    
    if (!prim.has_value()) return false;

    auto p = prim.value();
    if(!features.enableTextureMapping){
        hitInfo.material = *p.mat;
    }else{
        if(*p.mat->kdTexture){
            
        }
    }

    hitInfo.normal = std::visit(
        make_visitor(
        [&](const TrianglePrim& t) {
            auto v1 = t.v2->position - t.v1->position;
            auto v2 = t.v3->position - t.v1->position;
            return glm::normalize(glm::cross(v1, v2));
        },
        [&](const SpherePrim& s) {
            auto p = ray.origin + ray.direction * ray.t;
            return glm::normalize(p - *s.c);
        }
        ), p.p);

    return prim.has_value();
}