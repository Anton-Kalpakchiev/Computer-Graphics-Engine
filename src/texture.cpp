#include "texture.h"
#include <cmath>
#include <framework/image.h>

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

    int width = image.width - 1;
    int height = image.height - 1;
    float i = width * texCoord.x + 0.5;
    float j = height * texCoord.y + 0.5;

    float min_X = 0.0f;
    float max_X = 0.0f;

    float floor_X = floor(i);
    float diff_X = i - floor_X;
    if(diff_X > 0.5f){
        min_X = floor_X + 0.5;
        max_X = ceil(i) + 0.5;
    }else{
        min_X = floor_X - 0.5;
        max_X = floor_X + 0.5;
    }

    float min_Y = 0.0f;
    float max_Y = 0.0f;

    float floor_Y = floor(j);
    float diff_Y = j - floor_Y;
    if(diff_Y > 0.5f){
        min_Y = floor_Y + 0.5;
        max_Y = ceil(j) + 0.5;
    }else{
        min_Y = floor_Y - 0.5;
        max_Y = floor_Y + 0.5;
    }

    if(fabs(i - min_X) < fabs(max_X - i) && fabs(j - min_Y) < fabs(max_Y - j)){
        int pixel_X = int(min_X - 0.5);
        int pixel_Y = int(min_Y - 0.5);
        return image.pixels[pixel_Y * width + pixel_X];
    }else if(fabs(i - min_X) < fabs(max_X - i) && fabs(j - min_Y) > fabs(max_Y - j)){
        int pixel_X = int(min_X - 0.5);
        int pixel_Y = int(min_Y - 0.5);
        return image.pixels[pixel_Y * width + pixel_X];
    }

}