#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <iostream>
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif


glm::vec3 recursiveRayTrace(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth, int rayDepthInitial){
    HitInfo hitInfo;
    bvh.setRecursionLevel(rayDepthInitial - rayDepth);
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            if(!(reflection.direction == glm::vec3(0.0f) && reflection.origin == glm::vec3(0.0f) && reflection.t == 0.0f)){
                    if(rayDepth > 0){
                        Lo += recursiveRayTrace(scene, bvh, reflection, features, rayDepth - 1, rayDepthInitial);
                    }
            }
        }

        // Draw a white debug ray if the ray hits.
        if(features.enableShading){
            drawRay(ray, Lo);
        }else{
            drawRay(ray, glm::vec3(1.0f));
        }
        
        float hardShadowAverage = 0.0f;
        for(const auto& light : scene.lights){
            if (std::holds_alternative<PointLight>(light)) {
                    const PointLight& pointLight = std::get<PointLight>(light);
                    hardShadowAverage += testVisibilityLightSample(pointLight.position, glm::vec3 {1.f}, bvh, features, ray, hitInfo);
            } 
        }
        if(scene.lights.size() > 0){
            if(hardShadowAverage > 0.0f){
                Lo *= 1.0f;
            }else{
                Lo *= 0.0f;
            }
        }

        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    return recursiveRayTrace(scene, bvh, ray, features, rayDepth, rayDepth);
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features, 5));
        }
    }
}