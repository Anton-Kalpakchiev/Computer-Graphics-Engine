#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>

extern float bloomScalar;
extern float bloomThreshold;
extern int bloomDebugOption;

// Forward declarations.

struct Scene;
class Screen;
class Trackball;
class BvhInterface;
struct Features;

// Main rendering function.
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);

// Get the color of a ray.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth = 0);

glm::mat3 weightsGaussian(float sigma);