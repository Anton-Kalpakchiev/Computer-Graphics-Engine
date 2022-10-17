cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(ComputerGraphics C CXX)

option(USE_PREBUILT_INTERSECT "Enable using prebuilt intersection library" ON)

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/framework")
	# Create framework library and include CMake scripts (compiler warnings, sanitizers and static analyzers).
	add_subdirectory("framework")
	set(REFERENCE_MODE OFF)
else()
	# During development the framework lives in parent folder.
	add_subdirectory("../../../framework/" "${CMAKE_BINARY_DIR}/framework/")
	set(REFERENCE_MODE ON)
endif()

find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(FinalProjectLib
	"src/scene.cpp"
	"src/draw.cpp"
	"src/screen.cpp"
	"src/bounding_volume_hierarchy.cpp"
	"src/bvh_interface.cpp"
	"src/light.cpp"
	"src/config.cpp"
	"src/texture.cpp"
	"src/shading.cpp"
	"src/interpolate.cpp"
	"src/render.cpp"
)

if (REFERENCE_MODE)
	target_sources(FinalProjectLib PRIVATE
		"src/extra/motion_blur.cpp"
		"src/extra/multiple_rays.cpp"
		"src/extra/bloom_effect.cpp"
		"src/extra/glossy_reflection.cpp"
	)
endif()

target_include_directories(FinalProjectLib PUBLIC "src")
target_link_libraries(FinalProjectLib PUBLIC CGFramework OpenGL::GLU OpenMP::OpenMP_CXX)
target_compile_features(FinalProjectLib PUBLIC cxx_std_20)
enable_sanitizers(FinalProjectLib)
set_project_warnings(FinalProjectLib)

if (USE_PREBUILT_INTERSECT)
	add_library(Intersect2 STATIC IMPORTED)
	if (WIN32)
		set_target_properties(Intersect2 PROPERTIES IMPORTED_LOCATION_DEBUG "${CMAKE_SOURCE_DIR}/prebuilt/Intersect_x64_debug.lib")
		set_target_properties(Intersect2 PROPERTIES IMPORTED_LOCATION_RELEASE "${CMAKE_SOURCE_DIR}/prebuilt/Intersect_x64_release.lib")
	elseif (APPLE)
		set_target_properties(Intersect2 PROPERTIES IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/prebuilt/libIntersect_mac_arm64.a")
	else()
		set_target_properties(Intersect2 PROPERTIES IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/prebuilt/libIntersect_linux_x64.a")
	endif()
	target_link_libraries(FinalProjectLib PUBLIC Intersect2)
else()
	target_sources(FinalProjectLib PRIVATE "src/intersect.cpp")
endif()

target_compile_definitions(FinalProjectLib PUBLIC
	"-DDATA_DIR=\"${CMAKE_CURRENT_LIST_DIR}/data/\"")

add_executable(FinalProject	"src/main.cpp")
target_link_libraries(FinalProject PUBLIC FinalProjectLib)
target_compile_features(FinalProject PUBLIC cxx_std_20)
enable_sanitizers(FinalProject)
set_project_warnings(FinalProject)
