#pragma once
#include "scene.h"
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <utility>
#include <variant>
#include <vector>
#include "common.h"

struct CameraConfig {
    float fieldOfView = 50.0f; // in degrees
    float distanceFromLookAt = 3.0f;
    glm::vec3 lookAt = { 0.0f, 0.0f, 0.0f };
    glm::vec3 rotation = { 20.0f, 20.0f, 0.0f }; // in degrees
};

struct Config {
    Features features = {};

    bool cliRenderingEnabled = false;
    glm::ivec2 windowSize = { 800, 800 };
    std::filesystem::path dataPath = DATA_DIR;
    std::variant<SceneType, std::filesystem::path> scene = SceneType::SingleTriangle;
    std::filesystem::path outputDir = "";
    std::vector<CameraConfig> cameras;
    std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight>> lights;
};

std::ostream& operator<<(std::ostream& arg, const Config& config);

Config readConfigFile(const std::filesystem::path& config_path);

std::string serialize(const SceneType& sceneType);
std::optional<SceneType> deserialize(const std::string& lowered);