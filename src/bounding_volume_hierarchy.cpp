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
#include <stack>
#include <iostream>

glm::vec3 triangleCenter(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
    return a + b + c / 3.f;
}

std::optional<AxisAlignedBox> getBoundingBox(const std::span<Primitive>& primitives, size_t beg, size_t end, const std::span<Vertex>& vertices) {
    std::optional<AxisAlignedBox> res;
    
    for (auto it = primitives.begin() + beg; it != primitives.begin() + end; it++) {
        float xMin, xMax, yMin, yMax, zMin, zMax;

        if (it->p.index() == 0) { // check if it's a triangle
            auto t = std::get<Triangle>(it->p);
            xMin = std::min({ vertices[t.x].position.x, vertices[t.y].position.x, vertices[t.z].position.x });
            xMax = std::max({ vertices[t.x].position.x, vertices[t.y].position.x, vertices[t.z].position.x });
            yMin = std::min({ vertices[t.x].position.y, vertices[t.y].position.y, vertices[t.z].position.y });
            yMax = std::max({ vertices[t.x].position.y, vertices[t.y].position.y, vertices[t.z].position.y });
            zMin = std::min({ vertices[t.x].position.z, vertices[t.y].position.z, vertices[t.z].position.z });
            zMax = std::max({ vertices[t.x].position.z, vertices[t.y].position.z, vertices[t.z].position.z });
        } else if (it->p.index() == 1) { // check if it's a sphere
            auto s = std::get<Sphere>(it->p);
            xMin = s.center.x - s.radius;
            xMax = s.center.x + s.radius;
            yMin = s.center.y - s.radius;
            yMax = s.center.y + s.radius;
            zMin = s.center.z - s.radius;
            zMax = s.center.z + s.radius;
        }

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
    m_numLevels = std::max(m_numLevels, (int)depth);
    auto aabb = getBoundingBox(primitives, beg, end, vertices).value();
    if (beg + 1 == end) {
        nodes.push_back(Node { aabb, 0, 0, beg, end, depth });
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
    nodes.push_back(Node {aabb, left, right, beg, end, depth });
    return nodes.size() - 1;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    
    m_numLevels = 0;

    // Get all triangles of the scene into the BVH
    for (const auto& mesh : pScene->meshes) {

        // The vertex indices of the new triangles will be offset by the vertices already in the BVH
        const auto offset = vertices.size();
        std::copy(mesh.vertices.begin(), mesh.vertices.end(), std::back_inserter(vertices));
        std::transform(mesh.triangles.begin(), mesh.triangles.end(), std::back_inserter(primitives), [&](auto t) {
            auto p1 = mesh.vertices[t.x].position;
            auto p2 = mesh.vertices[t.y].position;
            auto p3 = mesh.vertices[t.z].position;
            return Primitive{ glm::uvec3(t.x + offset, t.y + offset, t.z + offset),
                              triangleCenter(p1, p2, p3) };
        });
    }
    // Get all spheres of the scene into the BVH
    std::transform(pScene->spheres.begin(), pScene->spheres.end(), std::back_inserter(primitives), [](auto s) {
        return Primitive{ Sphere(s), glm::vec3(s.center) };
    });

    // We have all the primitives and their centers in the primitves vector
    // Create the BVH itself
    root = createBVH(0, primitives.size(), 0, 0);

    m_numLeaves = primitives.size();
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return std::min(m_numLevels, 12);
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
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    for (const auto& node : nodes) {
        if (node.level == level) {
            drawAABB(node.aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
        }
    }
    // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
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
        bool hit = false;

        std::stack<Node> stack = std::stack<Node>();
        stack.push(nodes[root]);

        std::vector<float> t_vec = std::vector<float>();

        while(!stack.empty()){

            Node parent = stack.top();
            stack.pop();

            float rollBack = ray.t;

            if(parent.left == 0 && parent.right == 0){

                Primitive p = primitives[parent.beg];
                if(p.p.index() == 0){
                    auto t = std::get<Triangle>(p.p);
                    if(intersectRayWithTriangle(vertices[t.x].position, vertices[t.y].position, vertices[t.z].position, ray, hitInfo)){
                        t_vec.push_back(ray.t);
                        ray.t = rollBack;
                    }
                }else{
                    auto s = std::get<Sphere>(p.p);
                    if(intersectRayWithShape(s, ray, hitInfo)){
                        t_vec.push_back(ray.t);
                        ray.t = rollBack;
                    }
                }

            }else{

                Node left = nodes[parent.left];
                Node right = nodes[parent.right];

                float t_left = -1.0f;
                float t_right = -1.0f;

                if(intersectRayWithShape(left.aabb, ray)){
                    t_left = ray.t;
                    ray.t = rollBack;
                }

                if(intersectRayWithShape(right.aabb, ray)){
                    t_right = ray.t;
                    ray.t = rollBack;
                }

                if(t_left > 0.0f && t_right > 0.0f){
                    if(t_left < t_right){
                        stack.push(right);
                        stack.push(left);
                    }else{
                        stack.push(left);
                        stack.push(right);
                }
                }else if(t_left > 0.0f){
                    stack.push(left);
                }else if(t_right > 0.0f){
                    stack.push(right);
                }

            }
        }

        if(t_vec.size() > 0){
            hit = true;
            ray.t = *min_element(t_vec.begin(), t_vec.end());
        }

        return hit;
    }
}