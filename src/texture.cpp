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

    if(!features.extra.enableBilinearTextureFiltering){
        int i = std::max(texCoord.x * image.width, 0.0f);
        int j = std::max((1.0f - texCoord.y) * image.height, 0.0f);

        i = std::min(i, image.width - 1);
        j = std::min(j, image.height - 1);

        return image.pixels[std::floor(j) * image.width + std::floor(i)];
    }

    int yUp = int(std::clamp((float)ceil((1 - texCoord.y) * image.height), 0.0f, (float)image.height));
    int xUp = int(std::clamp((float)ceil(texCoord.x * image.width), 0.0f, (float)image.width));

    int yDown = int(std::clamp((float)floor((1 - texCoord.y) * image.height), 0.0f, (float)image.height));
    int xDown = int(std::clamp((float)floor(texCoord.x * image.width), 0.0f, (float)image.width));

    float x_remainder = fmod(texCoord.x * image.width, 1);
    float y_remainder = fmod(texCoord.y * image.height, 1);

    glm::vec3 xdyd = image.pixels[yDown * image.width + xDown];
    glm::vec3 xdyu = image.pixels[yUp * image.width + xDown];
    glm::vec3 xuyd = image.pixels[yDown * image.width + xUp];
    glm::vec3 xuyu = image.pixels[yUp * image.width + xUp];

    glm::vec3 y_1 = (1 - y_remainder) * xdyd + y_remainder * xuyd;
    glm::vec3 y_2 = (1 - y_remainder) * xdyu + y_remainder * xuyu;

    return (1 - x_remainder) * y_1 + x_remainder * y_2;

}