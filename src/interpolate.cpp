#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 a = v1 - v0;
    glm::vec3 b = v2 - v0;
    glm::vec3 c = p - v0;
    float den = a.x * b.y - b.x * a.y;
    glm::vec3 result = { 0.0f,
        0.0f,
        0.0f };
    result.y = (c.x * b.y - b.x * c.y) / den;
    result.z = (a.x * c.y - c.x * a.y) / den;
    result.x = 1 - result.y - result.z;
    return result;
    
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    return glm::normalize((n0 * barycentricCoord.x + n1 * barycentricCoord.y + n2 * barycentricCoord.z) / 3.0f);
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
// TODO: implement this function.
    return glm::vec2(0.0);
}
