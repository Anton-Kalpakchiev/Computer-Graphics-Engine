#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>
#include <vector>
#include <common.h>

extern float bloomScalar;
extern float bloomThreshold;
extern int bloomDebugOption;
extern int raysPerPixelSide;
extern int raysPerReflection;
extern float alphaModifier;

// Forward declarations.

struct Scene;
class Screen;
class Trackball;
class BvhInterface;
struct Features;

extern int raysPerPixelSide;
extern int samplesDoF;
extern float focusPlaneDistance;
extern float blurStrength;

// Main rendering function.
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);

// Get the color of a ray.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth = 0);

glm::mat3 weightsGaussian(float sigma);
std::vector<Ray> getRaySamples(const glm::vec2& pixelPos, const glm::vec2& pixelSize, const Trackball& camera, int n);

std::vector<Ray> getDOFRays(const glm::vec2& pixelPos, const Trackball& camera, float focalLength, float samplingRadius, int n);

Plane getPlane(const Trackball& camera, float dist);