#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>


const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{

    if (!features.enableShading)
        return glm::vec3 { 0.0f, 0.0f, 0.0f };

    glm::vec3 newNormal = glm::normalize(hitInfo.normal);
    glm::vec3 light = glm::normalize(lightPosition - (ray.direction * ray.t + ray.origin));

    float d = glm::dot(newNormal, light);
    if (d < 0) {
        d = 0;
    }

    glm::vec3 diffuse = hitInfo.material.kd * lightColor * d;//difuse lighting finished

    //newNormal stays the same
    //light stays the same
    glm::vec3 camera = glm::normalize(ray.direction);

    d = 0.0f;

    if (glm::dot(newNormal, light) > 0 && glm::dot(newNormal, camera) > 0) {
        glm::vec3 reflection = 2.0f * glm::dot(light, newNormal) * newNormal - light;
        d = std::pow(glm::dot(camera, reflection), hitInfo.material.shininess);
    }

    glm::vec3 specular = hitInfo.material.ks * lightColor * d;//specular done

    return diffuse + specular;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    // TODO: implement the reflection ray computation.
    return reflectionRay;
}