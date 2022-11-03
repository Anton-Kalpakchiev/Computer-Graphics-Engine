#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <vector>

class Screen {
public:
    Screen(const glm::ivec2& resolution, bool presentable = true);

    void clear(const glm::vec3& color);
    void setPixel(int x, int y, const glm::vec3& color);
    std::vector<glm::vec3> getTextureData();
    void writeBitmapToFile(const std::filesystem::path& filePath);
    void draw();

    [[nodiscard]] glm::ivec2 resolution() const;

    /// Calculates the index of a pixel in the `m_textureData` vector.
    /// Pixels are stored from bottom to top, left to right (this is to facilitate writing as a bmp).
    [[nodiscard]] int indexAt(int x, int y) const;
    [[nodiscard]] const std::vector<glm::vec3>& pixels() const;
    [[nodiscard]] std::vector<glm::vec3>& pixels();

private:
    bool m_presentable;
    glm::ivec2 m_resolution;
    std::vector<glm::vec3> m_textureData;
    uint32_t m_texture;
};
