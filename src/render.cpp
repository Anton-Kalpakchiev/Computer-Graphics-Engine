#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <iostream>
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

float bloomScalar = .3f;
float bloomThreshold = .4f;
int bloomDebugOption = 0;//afterPicture



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

//Adds bloom filted to the final image
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
    for (int y = 0; y < windowResolution.y - 1; y++) { // compute boxfilter with gaussian distribution per pixel and add bloom
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
            if (bloomDebugOption == 0) {//showcase final image
                screen.setPixel(x, y, newColor);
            } else if (bloomDebugOption == 1) { // showcase addition of bloom
                screen.setPixel(x, y, sum * scalar);
            } else {//showcase image before
                screen.setPixel(x, y, screenData[idx]);
            }
        }
    }
}


glm::mat3 weightsGaussian(float sigma) {
    float sum = 0.0f;
    glm::mat3 answer = {};
    for (int i = -1; i < 2; i++) {
        for (int k = -1; k < 2; k++) {
            float weight = exp(-(i * i + k * k) / (2 * sigma * sigma)) / (2 * 3.1415 * sigma * sigma);
            answer[i+ 1][k + 1] = weight;
            sum += weight;
        }
    }
    return answer / sum;
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
    if (features.extra.enableBloomEffect) { 
        renderBloomFilter(screen, features);
    }
}