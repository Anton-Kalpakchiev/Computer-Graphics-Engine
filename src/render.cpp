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

int raysPerPixelSide = 1;
float focusPlaneDistance = 1.f;
float blurStrength = .1f;

glm::vec3 recursiveRayTrace(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth, int rayDepthInitial){
    HitInfo hitInfo;
    bvh.setRecursionLevel(rayDepthInitial - rayDepth);
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        /*float hardShadowAverage = 0.0f;
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
        }*/

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

std::vector<Ray> getRaySamples(const glm::vec2& pixelPos, const glm::vec2& pixelSize, const Trackball& camera, int n) {
    std::vector<Ray> res;

    glm::vec2 pixelBox = pixelSize / float(n);
    std::default_random_engine gen(time(NULL));
    std::uniform_real_distribution<float> xJitter(0.f, pixelBox.x);
    std::uniform_real_distribution<float> yJitter(0.f, pixelBox.y);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            glm::vec2 newPos = glm::vec2(pixelPos.x + i * pixelBox.x, pixelPos.y + j * pixelBox.y);
            Ray r = camera.generateRay(newPos + glm::vec2(xJitter(gen), yJitter(gen)));
            res.emplace_back(r);
        }
    }
    return res;
}

Plane getPlane(const Trackball& camera, float dist) {
    auto planeNormal = glm::normalize(camera.lookAt() - camera.position());
    auto point = camera.position() + planeNormal * dist;
    auto D = glm::dot(point, point);
    return Plane(D, planeNormal);
}

glm::vec3 getIntersection(const Ray& ray, const Plane& plane) {
    const auto& [D, n] = plane;
    const auto& [o, d, _] = ray;
    float t = (D - glm::dot(n, o)) / glm::dot(n, d);
    return o + t * d;
}

std::vector<Ray> getDOFRays(const glm::vec2& pixelPos, const Trackball& camera, float focalLength, float jitter, int n) {
    auto focalPlane = getPlane(camera, focalLength);
    auto cameraPlane = getPlane(camera, 0.f);
    auto ray = camera.generateRay(pixelPos);

    // v1 and v2 are the basis the plane
    auto N = cameraPlane.normal;
    auto v1 = glm::normalize(glm::vec3 { -N.y, N.x, 0});
    if (N.x == 0.f && N.y == 0.f) v1 = glm::normalize(glm::vec3 { N.z, 0, -N.x });
    auto v2 = glm::normalize(glm::cross(N, v1));

    auto focalPoint = getIntersection(ray, focalPlane);

    float side = jitter / n;

    std::default_random_engine gen(time(NULL));
    std::uniform_real_distribution<float> jitterGen(0.f, side);
    auto point = ray.origin - jitter * v1 / 2.f - jitter * v2 / 2.f;

    std::vector<Ray> res;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            auto newOrigin = point + i * side * v1 + j * side * v2;
            newOrigin += jitterGen(gen) * v1 + jitterGen(gen);
            auto direction = focalPoint - newOrigin;
            res.emplace_back(Ray(newOrigin, direction));
        }
    }
    return res;
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

            auto colorSum = glm::vec3(0.f);
            size_t raysCast = 0;

            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const glm::vec2 pixelSize {
                1 / float(windowResolution.x) * 2.f,
                1 / float(windowResolution.y) * 2.f
            };

            if (features.extra.enableMultipleRaysPerPixel) {
                for (auto& ray : getRaySamples(normalizedPixelPos, pixelSize, camera, raysPerPixelSide)) {
                    colorSum += getFinalColor(scene, bvh, ray, features, 5);
                }
                raysCast += raysPerPixelSide * raysPerPixelSide;
            }

            if (features.extra.enableDepthOfField) {
                for (auto& ray : getDOFRays(normalizedPixelPos, camera, focusPlaneDistance, blurStrength, raysPerPixelSide)) {
                    colorSum += getFinalColor(scene, bvh, ray, features, 5);
                }
                raysCast += raysPerPixelSide * raysPerPixelSide;
            }

            if (!features.extra.enableMultipleRaysPerPixel) {
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                colorSum += getFinalColor(scene, bvh, cameraRay, features, 5);
                raysCast++;
            }

            glm::vec3 finalColor = colorSum / float(raysCast);
            screen.setPixel(x, y, finalColor);
        }
    }
}