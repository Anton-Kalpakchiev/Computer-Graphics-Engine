#include "image.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
DISABLE_WARNINGS_POP()
#include <cassert>
#include <exception>
#include <iostream>
#include <string>
#include <algorithm>

Image::Image(const std::filesystem::path& filePath)
{
	if (!std::filesystem::exists(filePath)) {
		std::cerr << "Texture file " << filePath << " does not exists!" << std::endl;
		throw std::exception();
	}

	const auto filePathStr = filePath.string(); // Create l-value so c_str() is safe.
	[[maybe_unused]] int numChannelsInSourceImage;
	stbi_uc* stbPixels = stbi_load(filePathStr.c_str(), &width, &height, &numChannelsInSourceImage, STBI_rgb);

	if (!stbPixels) {
		std::cerr << "Failed to read texture " << filePath << " using stb_image.h" << std::endl;
		throw std::exception();
	}

	constexpr size_t numChannels = 3; // STBI_rgb == 3 channels
	for (size_t i = 0; i < width * height * numChannels; i += numChannels) {
            pixels.emplace_back(stbPixels[i + 0] / 255.0f, stbPixels[i + 1] / 255.0f, stbPixels[i + 2] / 255.0f);
	}

	int k = int(std::floor((float)log2(pixels.size()) / 2));
	levels = std::vector<std::vector<glm::vec3> >();
	levels.push_back(pixels);

	for(int x = 0; x < k; x++){

		int w = (width / pow(2, x));
		std::vector<glm::vec3> vec = std::vector<glm::vec3>();
		std::copy(levels[x].begin(), levels[x].end(), std::back_inserter(vec));
		std::vector<glm::vec3> v = std::vector<glm::vec3>();
		for(int i = 0; i < w - 1; i += 2){
			for(int j = 0; j < w - 1; j += 2){
				glm::vec3 v1 = vec[j * w + i];
				glm::vec3 v2 = vec[(j + 1) * w + i];
				glm::vec3 v3 = vec[j * w + (i + 1)];
				glm::vec v4 = vec[(j+1) * w + (i + 1)];
				glm::vec3 avg = (v1 + v2 + v3 + v4) / 4.0f;
				v.push_back(avg);
			}
		}
		levels.push_back(v);
	}

	stbi_image_free(stbPixels);
}


void Image::setLOD(int lod){
	this->lod = lod;
}