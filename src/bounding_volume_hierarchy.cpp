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

std::optional<AxisAlignedBox> getBoundingBox(std::vector<Primitive>::const_iterator beg, std::vector<Primitive>::const_iterator end, Scene* scene) {
    std::optional<AxisAlignedBox> res;
    
    for (; beg != end; beg++) {
        float xMin, xMax, yMin, yMax, zMin, zMax;

        std::visit(make_visitor(
            [&](const TrianglePrim& t) {
                const auto& vertices = scene->meshes[t.meshIdx].vertices;
                const auto& p1 = vertices[t.v1].position;
                const auto& p2 = vertices[t.v2].position;
                const auto& p3 = vertices[t.v3].position;
                xMin = std::min({ p1.x, p2.x, p3.x });
                xMax = std::max({ p1.x, p2.x, p3.x });
                yMin = std::min({ p1.y, p2.y, p3.y });
                yMax = std::max({ p1.y, p2.y, p3.y });
                zMin = std::min({ p1.z, p2.z, p3.z });
                zMax = std::max({ p1.z, p2.z, p3.z });
            },
            [&](const SpherePrim& sp) {
                const auto& s = scene->spheres[sp.sphereIdx]; 
                xMin = s.center.x - s.radius;
                xMax = s.center.x + s.radius;
                yMin = s.center.y - s.radius;
                yMax = s.center.y + s.radius;
                zMin = s.center.z - s.radius;
                zMax = s.center.z + s.radius;
            }
        ), beg->p);

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


inline float boundingBoxSurfaceArea(const AxisAlignedBox& box) {
    auto lengths = box.upper - box.lower;
    return 2 * (lengths.x * lengths.y + lengths.y * lengths.z + lengths.z * lengths.x);
}

inline glm::vec3 triangleCenter(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
    return a + b + c / 3.f;
}

size_t splitStandard(std::vector<Primitive>& prims, size_t beg, size_t end, size_t depth, Scene*, std::vector<SAHCuts>&) {
    size_t mid = beg + (end - beg) / 2;
    std::nth_element(prims.begin() + beg, prims.begin() + mid, prims.begin() + end, comparators[depth % 3]);
    return mid;
}

float calculateSplitCost(std::vector<Primitive>& prims, size_t beg, size_t end, size_t split, Scene* scene) {
    auto boxLeft = getBoundingBox(prims.begin() + beg, prims.begin() + split, scene).value_or(AxisAlignedBox { glm::vec3(0.f), glm::vec3(0.f)});
    auto boxRight = getBoundingBox(prims.begin() + split, prims.begin() + end, scene).value_or(AxisAlignedBox { glm::vec3(0.f), glm::vec3(0.f)});

    auto areaLeft = boundingBoxSurfaceArea(boxLeft);
    auto areaRight = boundingBoxSurfaceArea(boxRight);
    
    return areaLeft * (split - beg) + areaRight * (end - split);
}

AxisAlignedBox getSplitPlane(glm::vec3 splitPoint, const AxisAlignedBox& container, size_t axis) {
    auto [low, high] = container;
    low[axis] = splitPoint[axis];
    high[axis] = splitPoint[axis];
    return AxisAlignedBox(low, high);
}

size_t splitSAHBinning(std::vector<Primitive>& prims, size_t beg, size_t end, size_t depth, Scene* scene, std::vector<SAHCuts>& debugCuts) {

    size_t skip = std::max(1UL, (unsigned long) (end - beg) / (unsigned long) (NUM_OF_BINS));

    auto parentBox = getBoundingBox(prims.begin() + beg, prims.begin() + end, scene).value_or(AxisAlignedBox(glm::vec3(0), glm::vec3(0)));
    size_t bestSplit;
    size_t bestAxis;
    float bestCost = std::numeric_limits<float>::max();

    auto cuts = SAHCuts{};

    for (size_t axis = 0; axis < 3; axis++) {
        std::sort(prims.begin() + beg, prims.begin() + end, comparators[axis]);
        for (size_t split = beg + skip; split < end; split += skip) {
            auto cost = calculateSplitCost(prims, beg, end, split, scene);
            if (cost < bestCost) {
                bestSplit = split;
                bestAxis = axis;
                bestCost = cost;
            }
        }
    }

    std::sort(prims.begin() + beg, prims.begin() + end, comparators[bestAxis]);

    return bestSplit;
}

size_t BoundingVolumeHierarchy::createBVH(size_t beg, size_t end, size_t depth) {
    m_numLevels = std::max(m_numLevels, (int)depth + 1);
    auto aabb = getBoundingBox(primitives.begin() + beg, primitives.begin() + end, m_pScene).value();
    if (depth + 1 > sahCutsPerLevel.size()) {
        sahCutsPerLevel.resize(depth + 1);
    }
    if (depth + 1 == MAX_DEPTH || beg + 1 == end) {
        nodes.push_back(Node { aabb, { true, depth, beg, end } });
        m_numLeaves++;
        return nodes.size() - 1;
    }
    auto mid = splitFunc(primitives, beg, end, depth, m_pScene, sahCutsPerLevel[depth]);
    auto left = createBVH(beg, mid, depth + 1);
    auto right = createBVH(mid, end, depth + 1);
    nodes.push_back(Node {aabb, { false, depth, beg, end, left, right } });
    return nodes.size() - 1;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{

    // Get all triangles of the scene into the BVH
    for (size_t meshIdx = 0; const auto& mesh : pScene->meshes) {

        std::transform(mesh.triangles.begin(), mesh.triangles.end(), std::back_inserter(primitives), [&](const auto& t) {
            const auto& p1 = mesh.vertices[t.x];
            const auto& p2 = mesh.vertices[t.y];
            const auto& p3 = mesh.vertices[t.z];
            return Primitive { TrianglePrim { meshIdx, t.x, t.y, t.z }, triangleCenter(p1.position, p2.position, p3.position) };
        });
        meshIdx++;
    }
    // Get all spheres of the scene into the BVH
    for (size_t sphereIdx = 0; const auto& s : pScene->spheres) {
        primitives.push_back(Primitive { SpherePrim { sphereIdx }, glm::vec3(s.center) });
        sphereIdx++;
    }

    m_numLeaves = 0;
    m_numLevels = 0;
    m_recursionLevel = 0;
    RECURSION_LEVEL = -1;

    if (features.extra.enableBvhSahBinning) {
        splitFunc = splitSAHBinning;
    } else {
        splitFunc = splitStandard;
    }

    // We have all the primitives and their centers in the primitves vector
    // Create the BVH itself
    root = createBVH(0, primitives.size(), 0);
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

void BoundingVolumeHierarchy::setRecursionLevel(int level, bool debug){
    if(!debug){
        m_recursionLevel = level;
    }else{
        RECURSION_LEVEL = level;
    }
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

void BoundingVolumeHierarchy::debugDrawSAHSplits(int level, int axis) {

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
                        const auto& v1 = m_pScene->meshes[t.meshIdx].vertices[t.v1];
                        const auto& v2 = m_pScene->meshes[t.meshIdx].vertices[t.v2];
                        const auto& v3 = m_pScene->meshes[t.meshIdx].vertices[t.v3];
                        debugDrawTriangle(v1, v2, v3);
                    },
                    [&](const SpherePrim& s) {
                        debugDrawSphere(m_pScene->spheres[s.sphereIdx]);
                    }
                ), p.p);
            }
            break;
        }
    }
}

std::optional<Primitive> getIntersecting(auto beg, auto end, Ray& ray, HitInfo& hitInfo, Scene* scene) {
    std::optional<Primitive> res;
    auto range = std::ranges::subrange(beg, end);
    for (const auto& prim : range) {
        auto hit =  std::visit(make_visitor(
            [&](const TrianglePrim& t) {
                const auto& v1 = scene->meshes[t.meshIdx].vertices[t.v1];
                const auto& v2 = scene->meshes[t.meshIdx].vertices[t.v2];
                const auto& v3 = scene->meshes[t.meshIdx].vertices[t.v3];
                return intersectRayWithTriangle(v1.position, v2.position, v3.position, ray, hitInfo);
            },
            [&](const SpherePrim& s) {
                return intersectRayWithShape(scene->spheres[s.sphereIdx], ray, hitInfo);
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
        
        prim = getIntersecting(primitives.begin(), primitives.end(), ray, hitInfo, m_pScene);
        
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
                auto maybePrim = getIntersecting(primitives.begin() + beg, primitives.begin() + end, ray, hitInfo, m_pScene);
                if (maybePrim.has_value()) {
                    prim = maybePrim.value();
                }

            } else {

                auto left = nodes[parent.data[4]];
                auto right = nodes[parent.data[5]];

                auto rollBack = ray.t;

                ray.t = std::numeric_limits<float>::max();
                bool leftBox = intersectRayWithShape(left.aabb, ray);
                if(leftBox){
                    if(m_recursionLevel == RECURSION_LEVEL){
                        drawAABB(left.aabb, DrawMode::Wireframe, glm::vec3(0.9f));
                    }
                }

                ray.t = std::numeric_limits<float>::max();
                bool rightBox = intersectRayWithShape(right.aabb, ray);
                if(rightBox){
                    if(m_recursionLevel == RECURSION_LEVEL){
                        drawAABB(right.aabb, DrawMode::Wireframe, glm::vec3(0.9f));
                    }
                }

                ray.t = rollBack;
                
                if (leftBox) stack.push(left);
                if (rightBox) stack.push(right);
                if(!leftBox && !rightBox){
                    if(m_recursionLevel == RECURSION_LEVEL){
                        drawAABB(parent.aabb, DrawMode::Wireframe, {0.9f, 0.0f, 0.0f});
                    }
                }
            }
        }
    }
    
    if (!prim.has_value()) return false;

    auto p = prim.value();

    hitInfo.normal = std::visit(make_visitor(
        [&](const TrianglePrim& t) {

            const auto& v1 = m_pScene->meshes[t.meshIdx].vertices[t.v1];
            const auto& v2 = m_pScene->meshes[t.meshIdx].vertices[t.v2];
            const auto& v3 = m_pScene->meshes[t.meshIdx].vertices[t.v3];

            if (m_recursionLevel == RECURSION_LEVEL && features.enableAccelStructure){
                debugDrawTriangle(v1, v2, v3);
            }

            if (features.enableNormalInterp) {
                glm::vec3 barCoords = computeBarycentricCoord(v1.position, v2.position, v3.position, ray.origin + ray.direction * ray.t);
                glm::vec3 interpolatedNormal = interpolateNormal(v1.normal, v2.normal, v3.normal, barCoords);
                if (glm::dot(interpolatedNormal, ray.direction) > 0) {
                    interpolatedNormal.x = -interpolatedNormal.x;
                    interpolatedNormal.y = -interpolatedNormal.y;
                    interpolatedNormal.z = -interpolatedNormal.z;
                }
                Ray toDraw = Ray(ray.origin + ray.direction * ray.t, interpolatedNormal, 1);
                drawRay(toDraw);
                drawRay(Ray(v1.position, v1.normal, 1));
                drawRay(Ray(v2.position, v2.normal, 1));
                drawRay(Ray(v3.position, v3.normal, 1));
                return interpolatedNormal;
            } else {
                auto u1 = v2.position - v1.position;
                auto u2 = v3.position - v1.position;
                return glm::normalize(glm::cross(u1, u2));
            }
        },
        [&](const SpherePrim& s) {
            auto p = ray.origin + ray.direction * ray.t;
            return glm::normalize(p - m_pScene->spheres[s.sphereIdx].center);
        }
    ), p.p);


    hitInfo.material = std::visit(make_visitor(
        [&](const TrianglePrim& t) {
            const auto& v1 = m_pScene->meshes[t.meshIdx].vertices[t.v1];
            const auto& v2 = m_pScene->meshes[t.meshIdx].vertices[t.v2];
            const auto& v3 = m_pScene->meshes[t.meshIdx].vertices[t.v3];
            auto material = m_pScene->meshes[t.meshIdx].material;
            if (!features.enableTextureMapping || !material.kdTexture) {
                return material;
            }
            const auto& texCoord = interpolateTexCoord(v1.texCoord, v2.texCoord, v3.texCoord, 
                                        computeBarycentricCoord(v1.position, v2.position, v3.position, ray.t * ray.direction + ray.origin));
            material.kd = acquireTexel(*material.kdTexture.get(), texCoord, features);
            return material;
        },
        [&](const SpherePrim& s) {
            return m_pScene->spheres[s.sphereIdx].material; 
        }    
    ), p.p);

    return true;
}