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

glm::vec3 triangleCenter(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
    return a + b + c / 3.f;
}

std::optional<AxisAlignedBox> getBoundingBox(std::span<Primitive> primitives, size_t beg, size_t end, std::span<Vertex> vertices) {
    std::optional<AxisAlignedBox> res;
    
    for (auto it = primitives.begin() + beg; it != primitives.begin() + end; it++) {
        float xMin, xMax, yMin, yMax, zMin, zMax;

        std::visit(make_visitor(
            [&](const Triangle& t) {
                const auto& v1 = vertices[t.vertexIdx.x].position;
                const auto& v2 = vertices[t.vertexIdx.y].position;
                const auto& v3 = vertices[t.vertexIdx.z].position;
                xMin = std::min({ v1.x, v2.x, v3.x });
                xMax = std::max({ v1.x, v2.x, v3.x });
                yMin = std::min({ v1.y, v2.y, v3.y });
                yMax = std::max({ v1.y, v2.y, v3.y });
                zMin = std::min({ v1.z, v2.z, v3.z });
                zMax = std::max({ v1.z, v2.z, v3.z });
            },
            [&](const Sphere& s) {
                xMin = s.center.x - s.radius;
                xMax = s.center.x + s.radius;
                yMin = s.center.y - s.radius;
                yMax = s.center.y + s.radius;
                zMin = s.center.z - s.radius;
                zMax = s.center.z + s.radius;
            }
        ), it->p);

        auto v = res.value_or(AxisAlignedBox { glm::vec3(xMin, yMin, zMin), glm::vec3(xMax, yMax, zMax) });
        v.lower.x = std::min(v.lower.x, xMin);
        v.lower.y = std::min(v.lower.y, yMin);
        v.lower.z = std::min(v.lower.z, zMin);
        v.upper.x = std::max(v.upper.x, xMax);
        v.upper.y = std::max(v.upper.y, yMax);
        v.upper.z = std::max(v.upper.z, zMax);
        res = std::make_optional(v);
    }
    return res;
}

size_t BoundingVolumeHierarchy::createBVH(size_t beg, size_t end, size_t splitBy, size_t depth) {
    m_numLevels = std::max(m_numLevels, (int)depth + 1);
    auto aabb = getBoundingBox(primitives, beg, end, vertices).value();
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
    // Get all triangles of the scene into the BVH
    for (size_t i = 0; i < pScene->meshes.size(); i++) {
        const auto& mesh = pScene->meshes[i];

        // The vertex indices of the new triangles will be offset by the vertices already in the BVH
        const auto offset = vertices.size();
        std::copy(mesh.vertices.begin(), mesh.vertices.end(), std::back_inserter(vertices));
        std::transform(mesh.triangles.begin(), mesh.triangles.end(), std::back_inserter(primitives), [&](auto t) {
            auto v1 = mesh.vertices[t.x].position;
            auto v2 = mesh.vertices[t.y].position;
            auto v3 = mesh.vertices[t.z].position;
            return Primitive{ Triangle { glm::uvec3(t.x + offset, t.y + offset, t.z + offset), i },
                              triangleCenter(v1, v2, v3) };
        });
    }
    // Get all spheres of the scene into the BVH
    std::transform(pScene->spheres.begin(), pScene->spheres.end(), std::back_inserter(primitives), [](auto s) {
        return Primitive{ Sphere(s), glm::vec3(s.center) };
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
        if (leafIdx < 0) {
            drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(.2f, 1.f, .2f), .5f);
            setColor({1.f, 1.f, 0.f});
            for (const auto& p : std::ranges::subrange(primitives.begin() + node.data[2], primitives.begin() + node.data[3])) {
                std::visit(make_visitor(
                    [&](const Triangle& t) {
                        drawTriangle(vertices[t.vertexIdx.x], vertices[t.vertexIdx.y], vertices[t.vertexIdx.z]);
                    },
                    [&](const Sphere& s) {
                        drawSphere(s);
                    }
                ), p.p);
            }
            break;
        }
    }
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}