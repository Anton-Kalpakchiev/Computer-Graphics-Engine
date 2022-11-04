#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <fmt/printf.h>
#include <framework/trackball.h>
#include <iostream>
#include <random>
#include <common.h>
#ifdef NDEBUG
#include <omp.h>
#endif

int raysPerPixelSide = 3;//implemented as a slider for multipleRaysPerPixel
int samplesDoF = 5;
float focusPlaneDistance = 3.f;
float blurStrength = .005f;

float bloomScalar = .3f;//implemented as a slider
float bloomThreshold = .4f;//implemented as a slider
int bloomDebugOption = 0; //implemented as a slider, when 0 (default), shows the resulting picture

int glossyReflectionsCap = 3; // cap of glossy recursive reflections
int raysPerReflection = 40;//implemented as a slider
float alphaModifier = 1.f;//implemented as a slider

glm::vec3 recursiveRayTrace(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth, int rayDepthInitial)
{
    HitInfo hitInfo;
    bvh.setRecursionLevel(rayDepthInitial - rayDepth);
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
        Ray reflection = computeReflectionRay(ray, hitInfo);

        if(features.extra.enableTransparency && !features.enableRecursive && !features.extra.enableGlossyReflection){
            if(hitInfo.material.transparency == 1.0f || rayDepth <= 0){
                return Lo;
            }
            Lo *= hitInfo.material.transparency;
            Ray t = Ray();
            t.origin = (0.00001f + ray.t) * ray.direction + ray.origin;
            t.direction = ray.direction;
            HitInfo h = HitInfo();
            h.normal = hitInfo.normal;
            bvh.intersect(t, h, features);
            drawRay(t, Lo);
            Lo += (1.0f - hitInfo.material.transparency) * recursiveRayTrace(scene, bvh, t, features, rayDepth - 1, rayDepthInitial);
        }


        if (!((reflection.direction == glm::vec3(0.0f) && reflection.origin == glm::vec3(0.0f) && reflection.t == 0.0f) || rayDepth < 1)) {

            if (features.enableRecursive) {
                glm::vec3 originalDirection = reflection.direction; // this is reflection of r
                if (features.extra.enableGlossyReflection && hitInfo.material.shininess != 0) {
                    glm::vec3 w = glm::normalize(originalDirection);
                    glm::vec3 t = w;
                    float min = t.x;
                    int minIdx = 0;
                    if (t.y < min) {
                        min = t.y;
                        minIdx = 1;
                    }
                    if (t.z < min) {
                        min = t.z;
                        minIdx = 2;
                    }
                    t[minIdx] = 1.0f;
                    glm::vec3 u = glm::cross(t, w) / glm::length(glm::cross(t, w));
                    glm::vec3 v = glm::cross(w, u);
                    float a = (1 / hitInfo.material.shininess) * alphaModifier;

                    //Glossy reflections visual debugger
                    glm::vec3 v0 = reflection.origin + (w + (u / 2.0f + v / 2.0f) * a) * .7f;
                    glm::vec3 v1 = reflection.origin + (w + (u / 2.0f + v / 2.0f - v) * a)  * .7f;
                    glm::vec3 v2 = reflection.origin + (w + (u / 2.0f + v / 2.0f - v - u) * a)  * .7f;
                    glm::vec3 v3 = reflection.origin + (w + (u / 2.0f + v / 2.0f - u) * a) * .7f;
                    glm::vec3 debugColor = glm::vec3(1.0f, 0.0f, 1.0f);
                    float debugTransparency = .5f;
                    drawPlane(v0, v1, v2, v3, debugColor, debugTransparency);

                    glm::vec3 totalColor = { .0f, .0f, .0f };
                    for (int i = 0; i < raysPerReflection; i++) {
                        float randOne = float(rand()) / RAND_MAX;
                        float randTwo = float(rand()) / RAND_MAX;
                        float weightU = -a / 2 + randOne * a;
                        float weightV = -a / 2 + randTwo * a;
                        glm::vec3 glossReflection = w + weightU * u + weightV * v;
                        glossReflection = glm::normalize(glossReflection);
                        if (glm::dot(hitInfo.normal, glossReflection) > 0) { // angle is less than 90 degrees
                            Ray glossRay = Ray(reflection.origin, glossReflection, std::numeric_limits<float>::max());
                            glm::vec3 color = recursiveRayTrace(scene, bvh, glossRay, features, glm::min(rayDepth - 1, glossyReflectionsCap), rayDepthInitial);
                            totalColor += color * hitInfo.material.ks;
                        }
                    }
                    totalColor /= raysPerReflection;
                    Lo += totalColor;
                } else {
                    Lo += recursiveRayTrace(scene, bvh, reflection, features, rayDepth - 1, rayDepthInitial);
                }
            }

        }

        if (features.enableRecursive && !features.extra.enableGlossyReflection) {
            Ray reflection;
            if(features.extra.enableTransparency && hitInfo.material.transparency != 1.0f){
                reflection = Ray();
                reflection.origin = (0.00001f + ray.t) * ray.direction + ray.origin;
                reflection.direction = ray.direction;
            }else{
                reflection = computeReflectionRay(ray, hitInfo);
            }

            if(!(reflection.direction == glm::vec3(0.0f) && reflection.origin == glm::vec3(0.0f) && reflection.t == 0.0f)){
                    if(rayDepth > 0){
                        Lo += recursiveRayTrace(scene, bvh, reflection, features, rayDepth - 1, rayDepthInitial);
                    }
            }

            if(hitInfo.material.transparency != 1){
                glm::vec3 vec = recursiveRayTrace(scene, bvh, reflection, features, rayDepth - 1, rayDepthInitial);
                glm::vec3 v = hitInfo.material.transparency * vec + (1 - hitInfo.material.transparency) * Lo;
                HitInfo h = HitInfo();
                h.normal = hitInfo.normal;
                bvh.intersect(reflection, h, features);
                drawRay(ray, v);
                return v;
            }

        }


        // Draw a white debug ray if the ray hits.
        if (features.enableShading) {
            drawRay(ray, Lo);
        } else {
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

// Adds bloom filted to the final image
void renderBloomFilter(Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    std::vector<glm::vec3> screenData = screen.getTextureData(); // originalData
    std::vector<glm::vec3> screenThreshold = screen.getTextureData();
    float threshold = bloomThreshold;

    for (int i = 0; i < screenThreshold.size(); i++) { // assert threshhold
        float brightness = 0.2126 * screenThreshold[i].x + 0.7152 * screenThreshold[i].y + 0.0722 * screenThreshold[i].z;
        if (brightness < bloomThreshold) {
            screenThreshold[i] = glm::vec3 { 0.0f, 0.0f, 0.0f };
        }
    }
    for (int y = 0; y < windowResolution.y - 1; y++) { // compute filter with gaussian distribution per pixel and add bloom
        for (int x = 0; x < windowResolution.x - 1; x++) {
            int idx = screen.indexAt(x, y);
            glm::vec3 sum = { 0.0f, 0.0f, 0.0f };
            glm::mat3 weightMatrixGaussian = weightsGaussian(1.0f);
            for (int k = -1; k < 2; k++) {
                for (int j = -1; j < 2; j++) {
                    if (!(x + k < 0 || x + k > windowResolution.x - 1 || y + j < 0 || y + j > windowResolution.y - 1)) {
                        int thisIndex = screen.indexAt(x + k, y + j);
                        float weight = weightMatrixGaussian[k + 1][j + 1];
                        sum += screenThreshold[thisIndex] * weight;
                    }
                }
            }
            float scalar = bloomScalar;
            glm::vec3 newColor = screenData[idx] + sum * scalar;
            if (bloomDebugOption == 0) { // showcase final image
                screen.setPixel(x, y, newColor);
            } else if (bloomDebugOption == 1) { // showcase addition of bloom
                screen.setPixel(x, y, sum * scalar);
            } else { // showcase image before
                screen.setPixel(x, y, screenData[idx]);
            }
        }
    }
}

glm::mat3 weightsGaussian(float sigma)
{
    float sum = 0.0f;
    glm::mat3 answer = {};
    for (int i = -1; i < 2; i++) {
        for (int k = -1; k < 2; k++) {
            float weight = exp(-(i * i + k * k) / (2 * sigma * sigma)) / (2 * 3.1415 * sigma * sigma);
            answer[i + 1][k + 1] = weight;
            sum += weight;
        }
    }
    return answer / sum;
}
std::vector<Ray> getRaySamples(const glm::vec2& pixelPos, const glm::vec2& pixelSize, const Trackball& camera, int n) {
    std::vector<Ray> res;

    glm::vec2 pixelBox = pixelSize / float(n);
    std::random_device rd;
    std::mt19937 gen(rd());
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
    /*auto point = camera.position() + planeNormal * dist;*/
    /*auto D = glm::sqrt(glm::dot(point, point));*/
    return Plane(dist - (glm::sqrt(glm::dot(camera.position(), camera.position()))), planeNormal);
}

glm::vec3 getIntersection(const Ray& ray, const Plane& plane) {
    const auto& [D, n] = plane;
    const auto& [o, d, _] = ray;
    float t = (D - glm::dot(n, o)) / glm::dot(n, d);
    return o + t * d;
}

std::vector<Ray> getDOFRays(const glm::vec2& pixelPos, const Trackball& camera, float focalLength, float samplingRadius, int n) {
    auto focalPlane = getPlane(camera, focalLength);
    auto cameraPlane = getPlane(camera, 0.f);
    auto ray = camera.generateRay(pixelPos);

    // v1 and v2 are the basis the plane
    auto N = cameraPlane.normal;
    auto v1 = glm::normalize(glm::vec3 { -N.y, N.x, 0});
    if (N.x == 0.f && N.y == 0.f) v1 = glm::normalize(glm::vec3 { N.z, 0, -N.x });
    auto v2 = glm::normalize(glm::cross(N, v1));

    auto focalPoint = getIntersection(ray, focalPlane);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distrTheta(0.f, glm::pi<float>() * 2);
    std::uniform_real_distribution<float> distrR(0.f, samplingRadius);

    std::vector<Ray> res;
    for (int i = 0; i < n; i++) {
        float r = glm::sqrt(distrR(gen));
        float theta = distrTheta(gen);
        auto newOrigin = ray.origin + r * glm::cos(theta) * v1 + r * glm::sin(theta) * v2;
        auto newDirection = focalPoint - newOrigin;
        res.emplace_back(Ray(newOrigin, newDirection));
    }
    //auto point = ray.origin;

    //for (int i = 0; i < n; i++) {
        //for (int j = 0; j < n; j++) {
            //auto newOrigin = point + i * side * v1 + j * side * v2;
            //newOrigin += jitterGen(gen) * v1 + jitterGen(gen) * v2;
            //auto direction = focalPoint - newOrigin;
            //res.emplace_back(Ray(newOrigin, direction));
        //}
    //}
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
            size_t weight = 0;

            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const glm::vec2 pixelSize {
                1 / float(windowResolution.x) * 2.f,
                1 / float(windowResolution.y) * 2.f
            };

            if (features.extra.enableMultipleRaysPerPixel) {
                auto color = glm::vec3(0.f);
                for (auto& ray : getRaySamples(normalizedPixelPos, pixelSize, camera, raysPerPixelSide)) {
                    color += getFinalColor(scene, bvh, ray, features, 5);
                }
                color /= raysPerPixelSide * raysPerPixelSide;
                colorSum += color;
                weight++;
            }

            if (features.extra.enableDepthOfField) {
                auto color = glm::vec3(0.f);
                for (auto& ray : getDOFRays(normalizedPixelPos, camera, focusPlaneDistance, blurStrength, samplesDoF)) {
                    color += getFinalColor(scene, bvh, ray, features, 5);
                }
                color /= samplesDoF;
                // if both multiple rays per pixel and depth of field are enabled, we want the depth of field to affect the image more
                colorSum += color * 3.f;
                weight += 3;
            }

            if (!features.extra.enableMultipleRaysPerPixel && !features.extra.enableDepthOfField) {
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                colorSum += getFinalColor(scene, bvh, cameraRay, features, 5);
                weight++;
            }

            glm::vec3 finalColor = colorSum / float(weight);
            screen.setPixel(x, y, finalColor);
        }
    }
    if (features.extra.enableBloomEffect) {
        renderBloomFilter(screen, features);
    }
}