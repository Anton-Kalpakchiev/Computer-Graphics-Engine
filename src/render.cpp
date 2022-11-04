#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <fmt/printf.h>
#include <framework/trackball.h>
#include <iostream>
#include <random>
#ifdef NDEBUG
#include <omp.h>
#endif

int raysPerPixelSide = 1;//implemented as a slider for multipleRaysPerPixel

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
std::vector<Ray> getRaySamples(glm::vec2 pixelPos, glm::vec2 pixelSize, const Trackball& camera, int n)
{
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

            if (!features.extra.enableMultipleRaysPerPixel) {
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                colorSum += getFinalColor(scene, bvh, cameraRay, features, 5);
                raysCast++;
            }

            glm::vec3 finalColor = colorSum / float(raysCast);
            screen.setPixel(x, y, finalColor);
        }
    }
    if (features.extra.enableBloomEffect) {
        renderBloomFilter(screen, features);
    }
}