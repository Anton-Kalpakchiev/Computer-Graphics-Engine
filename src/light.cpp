#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <fmt/printf.h>


// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, 
                                const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (!features.enableHardShadow) return 1.0;
    // normalize the ray direction, recalculate t
    ray.t *= glm::length(ray.direction);
    ray.direction = glm::normalize(ray.direction);
    
    // add an offset to the ray to prevent self intersections
    auto p = ray.origin + ray.direction * (ray.t - .001f);

    float len = glm::length(samplePos - p);
    Ray toLight = Ray {p, glm::normalize(samplePos - p), std::numeric_limits<float>::max()};
    bool hit = bvh.intersect(toLight, hitInfo, features);
    if (hit && toLight.t < len) {
        drawRay(toLight, glm::vec3(1.f, 0.f, 0.f));
        return 0.f;
    } else {
        toLight.t = len;
        drawRay(toLight, glm::vec3(1.f));
        return 1.f;
    }
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        glm::vec3 result = { 0.0f, 0.0f, 0.0f };

        // If shading is enabled, compute the contribution from all lights.
        for (const auto& light : scene.lights) {
                 if (std::holds_alternative<PointLight>(light)) {
                     const PointLight pointLight = std::get<PointLight>(light);
                     glm::vec3 color = computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                     result += color;
                 } else if (std::holds_alternative<SegmentLight>(light)) {
                     const SegmentLight segmentLight = std::get<SegmentLight>(light);
                     // Perform your calculations for a segment light.
                 } else if (std::holds_alternative<ParallelogramLight>(light)) {
                     const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                     // Perform your calculations for a parallelogram light.
                 }
             }

        return result;
        // TODO: replace this by your own implementation of shading
        

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
