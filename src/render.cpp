#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <iostream>
#include <fmt/printf.h>
#include <random>
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

int raysPerPixel = 1;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            if(!(reflection.direction == glm::vec3(0.0f) && reflection.origin == glm::vec3(0.0f) && reflection.t == 0.0f)){
                    if(rayDepth > 0){
                        Lo += getFinalColor(scene, bvh, reflection, features, rayDepth - 1);
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

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    std::default_random_engine gen;
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {

            glm::vec3 finalColor;
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };

            if (!features.extra.enableMultipleRaysPerPixel) {
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                finalColor = getFinalColor(scene, bvh, cameraRay, features, 5);
            } else {
                const float xOffset = 1 / float(windowResolution.x) * 2.f;
                const float yOffset = 1 / float(windowResolution.y) * 2.f;
                std::uniform_real_distribution<float> xDistr(0.f, xOffset);
                std::uniform_real_distribution<float> yDistr(0.f, yOffset);

                auto colorSum = glm::vec3(0.f);
                for (size_t i = 0; i < raysPerPixel; i++) {
                    const auto newNormalizePixelPos = normalizedPixelPos + glm::vec2(xDistr(gen), yDistr(gen));
                    const auto cameraRay = camera.generateRay(newNormalizePixelPos);
                    colorSum += getFinalColor(scene, bvh, cameraRay, features, 5);
                }
                finalColor = colorSum / float(raysPerPixel);
            }
            screen.setPixel(x, y, finalColor);
        }
    }
}