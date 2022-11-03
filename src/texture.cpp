#include "texture.h"
#include <cmath>
#include <framework/image.h>
#include <cmath>
#include <iostream>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    if(!features.enableTextureMapping){
        return image.pixels[0];
    }

    /*if(!features.extra.enableMipmapTextureFiltering){
        int i = std::max(texCoord.x * image.width, 0.0f);
        int j = std::max((1.0f - texCoord.y) * image.height, 0.0f);

        i = std::min(i, image.width - 1);
        j = std::min(j, image.height - 1);

        return image.pixels[std::floor(j) * image.width + std::floor(i)];
    }*/

    auto v = image.levels[image.lod];
    int w = pow(2, log2(v.size()) / 2.0f);

    float i = std::max(texCoord.x * w, 0.0f);
    float j = std::max((1.0f - texCoord.y) * w, 0.0f);

    i = std::min(i, (float)w - 1);
    j = std::min(j, (float)w - 1);

    return v[std::floor(j) * w + std::floor(i)];
}