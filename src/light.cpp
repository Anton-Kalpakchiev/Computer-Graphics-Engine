#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <fmt/printf.h>
#include <iostream>

int segmentLightSamples = 25;
int parallelogramLightDirectionSamples = 5;



// samples a segment light source using jittering
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color, float index, float sampleSize)//index is from 0 to sampleSize
{
    float r = ((float)rand() / RAND_MAX);//random float between 0 and 1
    float weight = (index + r) / sampleSize;

    position = (segmentLight.endpoint1 - segmentLight.endpoint0) * weight + segmentLight.endpoint0;
    color = weight * segmentLight.color1 + (1 - weight) * segmentLight.color0;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color, float horizontalIndex, float verticalIndex, float sampleSizeA, float sampleSizeB)
{
    float horRandom = ((float)rand() / RAND_MAX);
    float verRandom = ((float)rand() / RAND_MAX);

    float horWeight = (horizontalIndex + horRandom) / sampleSizeA;
    float verWeight = (verticalIndex + verRandom) / sampleSizeB;

    glm::vec3 horVector = horWeight * parallelogramLight.edge01;
    glm::vec3 verVector = verWeight * parallelogramLight.edge02;
    position = parallelogramLight.v0 + horVector + verVector;

    glm::vec3 bottomColor = horWeight * parallelogramLight.color1 + (1 - horWeight) * parallelogramLight.color0;
    glm::vec3 topColor = horWeight * parallelogramLight.color3 + (1 - horWeight) * parallelogramLight.color2;
    color = verWeight * topColor + (1 - verWeight) * bottomColor;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, 
                                const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (!features.enableHardShadow && !features.enableSoftShadow) return 1.0;
    // normalize the ray direction, recalculate t
    ray.t *= glm::length(ray.direction);
    ray.direction = glm::normalize(ray.direction);
    
    // add an offset to the ray to prevent self intersections
    auto p = ray.origin + ray.direction * (ray.t - .00001f);

    Ray toLight = Ray {p, samplePos - p, 1.f};
    bool hit = bvh.intersect(toLight, hitInfo, features);
    if (hit) {
        drawRay(toLight, glm::vec3(1.f, 0.f, 0.f));
        return 0.f;
    } else {
        drawRay(toLight, debugColor);
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
                     float visibility = 1.0f;
                     if (features.enableHardShadow) {
                         visibility = testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo);
                     }
                     result += color * visibility;
                 } else if (std::holds_alternative<SegmentLight>(light)) {
                     const SegmentLight segmentLight = std::get<SegmentLight>(light);
                     if (features.enableSoftShadow) {
                         glm::vec3 color = { 0.0f, 0.0f, 0.0f };
                         float sampleSize = segmentLightSamples;
                         for (int i = 0; i < sampleSize; i++) {
                             glm::vec3 colorOfLight = { 0.0f, 0.0f, 0.0f };
                             glm::vec3 positionOfLight = { 0.0f, 0.0f, 0.0f };
                             sampleSegmentLight(segmentLight, positionOfLight, colorOfLight, i, sampleSize);
                             float visibility = testVisibilityLightSample(positionOfLight, colorOfLight, bvh, features, ray, hitInfo);
                             
                             glm::vec3 thisColor = computeShading(positionOfLight, colorOfLight, features, ray, hitInfo);
                             color += thisColor * visibility;
                         }
                         result += color / sampleSize;
                     }
                 } else if (std::holds_alternative<ParallelogramLight>(light)) {
                     const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                     if (features.enableSoftShadow) {
                         glm::vec3 color = { 0.0f, 0.0f, 0.0f };
                         float sampleSizeA = parallelogramLightDirectionSamples;
                         float sampleSizeB = parallelogramLightDirectionSamples;
                         for (int i = 0; i < sampleSizeA; i++) {
                             for (int k = 0; k < sampleSizeB; k++) {
                                 glm::vec3 colorOfLight = { 0.0f, 0.0f, 0.0f };
                                 glm::vec3 positionOfLight = { 0.0f, 0.0f, 0.0f };
                                 sampleParallelogramLight(parallelogramLight, positionOfLight, colorOfLight, i, k, sampleSizeA, sampleSizeB);
                                 float visibility = testVisibilityLightSample(positionOfLight, colorOfLight, bvh, features, ray, hitInfo);
                                 glm::vec3 thisColor = computeShading(positionOfLight, colorOfLight, features, ray, hitInfo);
                                 color += thisColor * visibility;
                             }
                         }
                         result += color / (sampleSizeA * sampleSizeB);
                     }
                 }
             }
        return result;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}