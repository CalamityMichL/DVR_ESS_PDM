/* Copyright Michael Rauter
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 the "License";
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "SimpleITK.h""
#include "sitkImageFileReader.h""

#include <memory>
#include <math.h>
#include <vector>
#include <iostream>
#include <string>
#include <format>
#include <windows.h>
#include <functional>

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include "../../libs_extern/toml/toml.hpp"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "../common/ShaderSpirv.h"
#include "../common/ComputeShaderSpirv.h"
#include "../common/TextureHelper.h"
#include "../common/Camera.h"
#include "../LaunchComputeShaders.h"
#include "../RenderGeometry.h"
#include "../Volume.h"
#include "../ITF.h"

#define FILENAME_INI "./config.ini"
#define DEFAULT_FILENAME_INPUT_VOLUME "./volumes/tooth.nrrd"
#define DEFAULT_FILETYPE_INPUT_VOLUME "NrrdImageIO"
#define DEFAULT_FILENAME_ITF "./itfs/default.itf"
#define DEFAULT_BLOCKSIZE 4
#define DEFAULT_NUM_INTENSITY_PARTITIONS 32
#define DEFAULT_USE_SPECIAL_PARTITION_FOR_ZERO_INTENSITY false
#define DEFAULT_NUM_TIMING_ITERATIONS 1

#define DEFAULT_BACKGROUND_COLOR_R 0.0f
#define DEFAULT_BACKGROUND_COLOR_G 0.0f
#define DEFAULT_BACKGROUND_COLOR_B 0.0f

#define PI 3.1415927f

std::string inputImageFileName = DEFAULT_FILENAME_INPUT_VOLUME;
std::string fileformat = DEFAULT_FILETYPE_INPUT_VOLUME;
std::string filename_itf = DEFAULT_FILENAME_ITF;
std::string filename_itf2 = DEFAULT_FILENAME_ITF;
bool force_volume_datatype_unorm8 = false;
bool precompute_gradient = true;
bool useFullScreen = false;
int enableTricubicFiltering = false;
int reductionBlockSize = DEFAULT_BLOCKSIZE;
unsigned int numIntensityPartitions = DEFAULT_NUM_INTENSITY_PARTITIONS;
int numTimingIterations = DEFAULT_NUM_TIMING_ITERATIONS;
bool useSpecialPartitionForZeroIntensity = DEFAULT_USE_SPECIAL_PARTITION_FOR_ZERO_INTENSITY;
bool useDemoAndTimingsMode = false;
double timingPeriodInSec = 10.0;
float rotationsPerTimingPeriod = 2.0;
float demoAndTimingsModeCameraDistance = 0.75f;
float demoAndTimingsModeCameraHeight = 0.0f;
float demoAndTimingsModeCameraPitch = 0.0f;
bool disablePartitionedDistanceMapOptimization;
float itfBlendTime = 3.0f;

int itfInterpolate = 0; // 0...don't interpolate, 1...interpolate to ITF1, 2...interpolate to ITF2
float t_itfInterpolation = 0.0f;

float initialVolumeRotationX = 0.0f;
float initialVolumeRotationY = 0.0f;
float initialVolumeRotationZ = 0.0f;

float background_color_r = DEFAULT_BACKGROUND_COLOR_R;
float background_color_g = DEFAULT_BACKGROUND_COLOR_G;
float background_color_b = DEFAULT_BACKGROUND_COLOR_B;

double PCFreq = 0.0;
__int64 CounterStart = 0;

void StartCounter()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		std::cout << "QueryPerformanceFrequency failed!\n";

	PCFreq = double(li.QuadPart) / 1000.0;

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}
double GetCounter()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}

template<class ReturnType, class... Ts>
ReturnType DoFunctionTimings(void* function, std::string text, int numTimingIterations, const Ts&... args) {
	reinterpret_cast<ReturnType(*)(Ts...)>(function)(args...);
	double summedTime = 0.0;
	for (int i = 0; i < numTimingIterations; i++)
	{
		glFlush();
		glFinish();
		StartCounter();
		reinterpret_cast<ReturnType(*)(Ts...)>(function)(args...);
		glFlush();
		glFinish();
		double elapsed_time = GetCounter();
		summedTime += elapsed_time;		
	}
	std::cout << text << " mean computation took (average of " << numTimingIterations << " iterations) " << summedTime / numTimingIterations << " msecs\n";
}

using namespace std;
namespace sitk = itk::simple;

Camera cam;

bool moveLeft, moveRight, moveForward, moveBackward = (false, false, false, false);
double lastTimeMovement, deltaTime = 0.0;
bool first_mouse = true;
double lastX, lastY;
bool first_movement = true;

int window_width = 1280;
int window_height = 720;
float fov = 60;

glm::mat4 projection;

struct ShaderParams {
	int ditheringMode = 2; // 0... dithering off, 1... procedurally generated, 2... noise texture
	int renderContentType = 0; // 0... volume rendering, 1... volume bounding box front faces, 2... volume bounding box back faces, 3... acceleration structure, 4... voxel gradient (only visualized, if lighting is incorporated into rendering)
	int essType = 2; // 0... no ess, 1... block ess, 2... distance map based ess
	int earlyRayTermination = 1;
	int gradientType = 0; // 0... on-the-fly, 1... use precomputed gradient
	int reuseNormals = 1;
	float samplingFactor = 1.0f;
	float gamma = 1.0f;
	float lightshadingFactor = 0.0f;
	float lightDirComponentX = -0.5f;
	float lightDirComponentY = 1.0f;
	float lightDirComponentZ = 0.25f;
	float lightAmbientFactor = 0.2f;
	float lightDiffuseFactor = 1.0f;
	float lightSpecularFactor = 1.0f;
	float lightSpecularShininess = 16.0f;
};
ShaderParams shaderParams;

GLuint textureNoise;
sitk::Image imgNoise;

void UpdateNoiseImage(sitk::Image &imgNoise)
{
	uint8_t* bufferNoise = imgNoise.GetBufferAsUInt8();
	for (unsigned int i = 0; i < imgNoise.GetWidth() * imgNoise.GetHeight(); i++)
		bufferNoise[i] = rand() % 256;
}

void WindowResizeCallback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	projection = glm::perspective((float)glm::radians(fov), (float)width / (float)height, 0.3f, 100.0f);
	window_width = width;
	window_height = height;

	sitk::Image imgNoise = sitk::Image(width, height, sitk::sitkUInt8);
	UpdateNoiseImage(imgNoise);
	glBindTexture(GL_TEXTURE_2D, textureNoise);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, imgNoise.GetBufferAsUInt8());

	first_mouse = true;
}

void KeyInputCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, 1);

	if (!useDemoAndTimingsMode)
	{
		if (key == GLFW_KEY_W && action == GLFW_PRESS)
			moveForward = true;
		else if (key == GLFW_KEY_W && action == GLFW_RELEASE)
			moveForward = false;
		if (key == GLFW_KEY_S && action == GLFW_PRESS)
			moveBackward = true;
		else if (key == GLFW_KEY_S && action == GLFW_RELEASE)
			moveBackward = false;
		if (key == GLFW_KEY_A && action == GLFW_PRESS)
			moveLeft = true;
		else if (key == GLFW_KEY_A && action == GLFW_RELEASE)
			moveLeft = false;
		if (key == GLFW_KEY_D && action == GLFW_PRESS)
			moveRight = true;
		else if (key == GLFW_KEY_D && action == GLFW_RELEASE)
			moveRight = false;
	}

	if (key == GLFW_KEY_1 && action == GLFW_PRESS)
		shaderParams.renderContentType = (shaderParams.renderContentType + 1) % 5;
	if (key == GLFW_KEY_2 && action == GLFW_PRESS)
		shaderParams.ditheringMode = (shaderParams.ditheringMode + 1) % 3;
	if (key == GLFW_KEY_3 && action == GLFW_PRESS)
		shaderParams.earlyRayTermination = (shaderParams.earlyRayTermination + 1) % 2;
	if (key == GLFW_KEY_4 && action == GLFW_PRESS)
		shaderParams.essType = (shaderParams.essType + 1) % 3;
	if (key == GLFW_KEY_5 && action == GLFW_PRESS)
	{
		shaderParams.lightshadingFactor = (float)((int)(shaderParams.lightshadingFactor * 10.0f) + 1) / 10.0f;
		if (shaderParams.lightshadingFactor > 1.0f)
			shaderParams.lightshadingFactor = 0.0f;
	}
	if (key == GLFW_KEY_6 && action == GLFW_PRESS)
		shaderParams.gradientType = (shaderParams.gradientType + 1) % 2;
	if (key == GLFW_KEY_7 && action == GLFW_PRESS)
		enableTricubicFiltering = (enableTricubicFiltering + 1) % 2;
	if (key == GLFW_KEY_8 && action == GLFW_PRESS)
		shaderParams.reuseNormals = (shaderParams.reuseNormals + 1) % 2;
	if (key == GLFW_KEY_9 && action == GLFW_PRESS)
	{
		shaderParams.lightDiffuseFactor = shaderParams.lightDiffuseFactor + 1.0;
		if (shaderParams.lightDiffuseFactor > 10.0f)
			shaderParams.lightDiffuseFactor = 1.0f;
	}
	if (key == GLFW_KEY_0 && action == GLFW_PRESS)
	{
		shaderParams.lightSpecularFactor = shaderParams.lightSpecularFactor + 1.0;
		if (shaderParams.lightSpecularFactor > 10.0f)
			shaderParams.lightSpecularFactor = 1.0f;
	}
	if (key == GLFW_KEY_F1 && action == GLFW_PRESS)
	{
		if ((itfInterpolate == 0) && (t_itfInterpolation == 1.0f))
		{
			itfInterpolate = 1;
			t_itfInterpolation = 1.0f;
		}
	}
	if (key == GLFW_KEY_F2 && action == GLFW_PRESS)
	{
		if ((itfInterpolate == 0) && (t_itfInterpolation == 0.0f))
		{
			itfInterpolate = 2;
			t_itfInterpolation = 0.0f;
		}
	}
}

void MouseLookCallback(GLFWwindow* window, double xpos, double ypos) {

	if (first_mouse) {
		lastX = xpos;
		lastY = ypos;
		first_mouse = false;
	}

	double xoffset = xpos - lastX;
	double yoffset = lastY - ypos;

	lastX = xpos;
	lastY = ypos;

	cam.ProcessMouseMovement(xoffset, yoffset, true);
}

void DoMovement() {
	if (first_movement) {
		lastTimeMovement = glfwGetTime();
		first_movement = false;
		return;
	}

	if (moveLeft)
		cam.ProcessKeyboard(Camera::Camera_Movement::LEFT, (float)deltaTime * 1.0f);
	if (moveRight)
		cam.ProcessKeyboard(Camera::Camera_Movement::RIGHT, (float)deltaTime * 1.0f);
	if (moveForward)
		cam.ProcessKeyboard(Camera::Camera_Movement::FORWARD, (float)deltaTime * 1.0f);
	if (moveBackward)
		cam.ProcessKeyboard(Camera::Camera_Movement::BACKWARD, (float)deltaTime * 1.0f);
	
	deltaTime = glfwGetTime() - lastTimeMovement;
	lastTimeMovement = glfwGetTime();
}

int main(int argc, char* argv[])
{
	std::string filenameIniFile = std::string(FILENAME_INI);
	if (argc > 1)
		filenameIniFile = std::string(argv[1]);

	toml::table tbl;
	try
	{
		tbl = toml::parse_file(filenameIniFile);
		if (auto optionalValue = tbl["settings"]["filename_input_volume"].value<string>())
			inputImageFileName = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["filetype_input_volume"].value<string>())
			fileformat = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["filename_itf"].value<string>())
			filename_itf = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["filename_itf2"].value<string>())
			filename_itf2 = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["initialVolumeRotation"][0].value<float>())
			initialVolumeRotationX = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["initialVolumeRotation"][1].value<float>())
			initialVolumeRotationY = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["initialVolumeRotation"][2].value<float>())
			initialVolumeRotationZ = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["force_volume_datatype_unorm8"].value<bool>())
			force_volume_datatype_unorm8 = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["precompute_gradient"].value<bool>())
			precompute_gradient = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["viewport_width"].value<int>())
			window_width = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["viewport_height"].value<int>())
			window_height = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["use_fullscreen"].value<bool>())
			useFullScreen = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["enableTricubicFiltering"].value<int>())
			enableTricubicFiltering = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["blocksize"].value<int>())
			reductionBlockSize = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["disable_partitioned_distancemap_optimization"].value<bool>())
			disablePartitionedDistanceMapOptimization = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["num_intensity_partitions"].value<unsigned int>())
			numIntensityPartitions = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["use_special_partition_for_zero_intensity"].value<bool>())
			useSpecialPartitionForZeroIntensity = optionalValue.value();	
		if (auto optionalValue = tbl["settings"]["use_demo_and_timings_mode"].value<bool>())
			useDemoAndTimingsMode = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["timing_period_in_sec"].value<double>())
			timingPeriodInSec = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["rotations_per_timing_period"].value<float>())
			rotationsPerTimingPeriod = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["demo_and_timings_mode_camera_distance"].value<float>())
			demoAndTimingsModeCameraDistance = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["demo_and_timings_mode_camera_height"].value<float>())
			demoAndTimingsModeCameraHeight = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["demo_and_timings_mode_camera_pitch"].value<float>())
			demoAndTimingsModeCameraPitch = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["num_timing_iterations"].value<int>())
			numTimingIterations = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["itf_blend_time"].value<float>())
			itfBlendTime = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["background_color"][0].value<float>())
			background_color_r = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["background_color"][1].value<float>())
			background_color_g = optionalValue.value();
		if (auto optionalValue = tbl["settings"]["background_color"][2].value<float>())
			background_color_b = optionalValue.value();
	}
	catch (const toml::parse_error& err)
	{
		std::cerr << "Parsing of ini-file failed:\n" << err << "\n";
		return 1;
	}

	try
	{
		lastX, lastY = window_width / 2, window_height / 2;

		if (!glfwInit())
			return(-1);

		GLFWmonitor* glfwMonitor = NULL;
		if (useFullScreen)
			glfwMonitor = glfwGetPrimaryMonitor();
		GLFWwindow* window = glfwCreateWindow(window_width, window_height, "Volume Rendering with C++", glfwMonitor, NULL);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);

		glfwMakeContextCurrent(window);

		if (!gladLoadGL(glfwGetProcAddress))
		{
			throw(std::string("Failed to initialize GLAD"));
		}

		// load shader render_volume program
		ShaderSpirv shader_render_volume("shaders/vs_render_volume.spv", "shaders/fs_render_volume.spv");
		// load compute shaders
		ComputeShaderSpirv computeshader_convert_volume_to_unorm8_from_int16 = ComputeShaderSpirv("shaders/cs_convert_volume_to_unorm8_from_dtype_int16.spv");
		ComputeShaderSpirv computeshader_convert_volume_to_unorm16_from_int16 = ComputeShaderSpirv("shaders/cs_convert_volume_to_unorm16_from_dtype_int16.spv");
		ComputeShaderSpirv computeshader_convert_volume_to_unorm8_from_uint16 = ComputeShaderSpirv("shaders/cs_convert_volume_to_unorm8_from_dtype_uint16.spv");
		ComputeShaderSpirv computeshader_convert_volume_to_unorm16_from_uint16 = ComputeShaderSpirv("shaders/cs_convert_volume_to_unorm16_from_dtype_uint16.spv");
		ComputeShaderSpirv computeshader_convert_volume_to_unorm8_from_uint8 = ComputeShaderSpirv("shaders/cs_convert_volume_to_unorm8_from_dtype_uint8.spv");
		ComputeShaderSpirv computeshader_convert_volume_to_unorm16_from_uint8 = ComputeShaderSpirv("shaders/cs_convert_volume_to_unorm16_from_dtype_uint8.spv");
		ComputeShaderSpirv computeshader_compute_gradient_map_unorm8 = ComputeShaderSpirv("shaders/cs_compute_gradient_map_dtype_unorm8.spv");
		ComputeShaderSpirv computeshader_compute_gradient_map_unorm16 = ComputeShaderSpirv("shaders/cs_compute_gradient_map_dtype_unorm16.spv");
		ComputeShaderSpirv computeshader_compute_occupancy_map_unorm8 = ComputeShaderSpirv("shaders/cs_compute_occupancy_map_dtype_unorm8.spv");
		ComputeShaderSpirv computeshader_compute_occupancy_map_unorm16 = ComputeShaderSpirv("shaders/cs_compute_occupancy_map_dtype_unorm16.spv");
		ComputeShaderSpirv computeshader_compute_POM_unorm8 = ComputeShaderSpirv("shaders/cs_compute_POM_dtype_unorm8.spv");
		ComputeShaderSpirv computeshader_compute_POM_unorm16 = ComputeShaderSpirv("shaders/cs_compute_POM_dtype_unorm16.spv");
		ComputeShaderSpirv computeshader_compute_distance_map = ComputeShaderSpirv("shaders/cs_compute_distance_map.spv");
		ComputeShaderSpirv computeshader_combine_PDMs = ComputeShaderSpirv("shaders/cs_combine_PDMs.spv");
		typedef ComputeShaderSpirv ComputeShaderType;
		ComputeShaderSpirv* computeshader_compute_gradient_map = NULL;
		ComputeShaderSpirv* computeshader_compute_occupancy_map = NULL;
		ComputeShaderSpirv* computeshader_compute_POM = NULL;

		RenderGeometry renderGeometry;

		glClearColor(background_color_r, background_color_g, background_color_b, 1.0f);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glEnable(GL_DEPTH_CLAMP);
		glEnable(GL_CULL_FACE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// disable vsync
		glfwSwapInterval(0);

		if (!useDemoAndTimingsMode)
		{
			glfwSetCursorPosCallback(window, MouseLookCallback);
		}
		glfwSetKeyCallback(window, KeyInputCallback);

		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

		std::vector<std::string> vecInputFileNames{ inputImageFileName };
		Volume volume(vecInputFileNames);

		vector<unsigned int> volumeSize = volume.GetVolumeSize();
		vector<double> voxelSize = volume.GetVoxelSize();
		vector<double> direction = volume.GetDirection();
		voxelSize[0] *= 0.001;
		voxelSize[1] *= 0.001;
		voxelSize[2] *= 0.001;

		itk::simple::PixelIDValueEnum pixelID = volume.GetPixelID();
		VoxelOutputDataType voxelOutputDataType;

		float minIntensityValue = (float)volume.minIntensityValue;
		float maxIntensityValue = (float)volume.maxIntensityValue;

		void* data = volume.GetBuffer();

		ITF itf;
		itf.LoadITF(filename_itf);

		ITF itf2;
		itf2.LoadITF(filename_itf2);

		ITF itfInterpolated;
		itfInterpolated.SetFromInterpolatedITFs(itf, itf2, t_itfInterpolation);

		uint32_t* itfData = itfInterpolated.GetData();

		sitk::Image imgNoise = sitk::Image(window_width, window_height, sitk::sitkUInt8);
		UpdateNoiseImage(imgNoise);

		int sizeOccupancyMap[3] = { 2 * ((int)(ceil((float)volumeSize[0] / (float)reductionBlockSize) + 1) / 2), 2 * ((int)(ceil((float)volumeSize[1] / (float)reductionBlockSize) + 1) / 2), 2 * ((int)(ceil((float)volumeSize[2] / (float)reductionBlockSize) + 1) / 2) };


		GLuint textureVolume;
		GLuint textureTmp;
		switch (pixelID)
		{
		case itk::simple::PixelIDValueEnum::sitkUInt16:
			if (!TextureHelper::CreateTexture3D(GL_R16UI, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED_INTEGER, GL_UNSIGNED_SHORT, (void*)data, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureTmp))
				throw std::runtime_error("error creating temporal volume texture");

			if (force_volume_datatype_unorm8)
			{
				if (!TextureHelper::CreateTexture3D(GL_R8, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED, GL_UNSIGNED_BYTE, NULL, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureVolume))
					throw std::runtime_error("error creating volume texture");

				ComputeConvertVolumeToSNORM(computeshader_convert_volume_to_unorm8_from_uint16, textureTmp, VoxelDataType::DTYPE_UINT16, VoxelOutputDataType::DTYPE_UNORM8, minIntensityValue, maxIntensityValue, textureVolume, volumeSize);
			}
			else
			{
				if (!TextureHelper::CreateTexture3D(GL_R16, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED, GL_UNSIGNED_SHORT, NULL, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureVolume))
					throw std::runtime_error("error creating volume texture");

				ComputeConvertVolumeToSNORM(computeshader_convert_volume_to_unorm16_from_uint16, textureTmp, VoxelDataType::DTYPE_UINT16, VoxelOutputDataType::DTYPE_UNORM16, minIntensityValue, maxIntensityValue, textureVolume, volumeSize);
			}

			TextureHelper::DeleteTexture(textureTmp);

			break;
		case itk::simple::PixelIDValueEnum::sitkInt16:
			if (!TextureHelper::CreateTexture3D(GL_R16I, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED_INTEGER, GL_SHORT, (void*)data, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureTmp))
				throw std::runtime_error("error creating temporal volume texture");

			if (force_volume_datatype_unorm8)
			{
				if (!TextureHelper::CreateTexture3D(GL_R8, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED, GL_UNSIGNED_BYTE, NULL, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureVolume))
					throw std::runtime_error("error creating volume texture");

				ComputeConvertVolumeToSNORM(computeshader_convert_volume_to_unorm8_from_int16, textureTmp, VoxelDataType::DTYPE_INT16, VoxelOutputDataType::DTYPE_UNORM8, minIntensityValue, maxIntensityValue, textureVolume, volumeSize);
			}
			else
			{
				if (!TextureHelper::CreateTexture3D(GL_R16, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED, GL_UNSIGNED_SHORT, NULL, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureVolume))
					throw std::runtime_error("error creating volume texture");

				ComputeConvertVolumeToSNORM(computeshader_convert_volume_to_unorm16_from_int16, textureTmp, VoxelDataType::DTYPE_INT16, VoxelOutputDataType::DTYPE_UNORM16, minIntensityValue, maxIntensityValue, textureVolume, volumeSize);
			}
			TextureHelper::DeleteTexture(textureTmp);

			break;
		case itk::simple::PixelIDValueEnum::sitkUInt8:
			if (!TextureHelper::CreateTexture3D(GL_R8UI, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED_INTEGER, GL_UNSIGNED_BYTE, (void*)data, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureTmp))
				throw std::runtime_error("error creating temporal volume texture");

			if (force_volume_datatype_unorm8)
			{
				if (!TextureHelper::CreateTexture3D(GL_R8, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED, GL_UNSIGNED_BYTE, NULL, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureVolume))
					throw std::runtime_error("error creating volume texture");

				ComputeConvertVolumeToSNORM(computeshader_convert_volume_to_unorm8_from_uint8, textureTmp, VoxelDataType::DTYPE_UINT8, VoxelOutputDataType::DTYPE_UNORM8, minIntensityValue, maxIntensityValue, textureVolume, volumeSize);
			}
			else
			{
				if (!TextureHelper::CreateTexture3D(GL_R16, volumeSize[0], volumeSize[1], volumeSize[2], GL_RED, GL_UNSIGNED_SHORT, NULL, GL_LINEAR, GL_CLAMP_TO_BORDER, &textureVolume))
					throw std::runtime_error("error creating volume texture");

				ComputeConvertVolumeToSNORM(computeshader_convert_volume_to_unorm16_from_uint8, textureTmp, VoxelDataType::DTYPE_UINT8, VoxelOutputDataType::DTYPE_UNORM16, minIntensityValue, maxIntensityValue, textureVolume, volumeSize);
			}
			TextureHelper::DeleteTexture(textureTmp);

			break;
		}

		if (force_volume_datatype_unorm8)
		{
			voxelOutputDataType = VoxelOutputDataType::DTYPE_UNORM8;
			computeshader_compute_gradient_map = &computeshader_compute_gradient_map_unorm8;
			computeshader_compute_occupancy_map = &computeshader_compute_occupancy_map_unorm8;
			computeshader_compute_POM = &computeshader_compute_POM_unorm8;
		}
		else
		{
			voxelOutputDataType = VoxelOutputDataType::DTYPE_UNORM16;
			computeshader_compute_gradient_map = &computeshader_compute_gradient_map_unorm16;
			computeshader_compute_occupancy_map = &computeshader_compute_occupancy_map_unorm16;
			computeshader_compute_POM = &computeshader_compute_POM_unorm16;
		}

		GLuint textureGradient;
		if (precompute_gradient)
		{
			if (!TextureHelper::CreateTexture3D(GL_RGBA8_SNORM, volumeSize[0], volumeSize[1], volumeSize[2], GL_RGBA, GL_BYTE, NULL, GL_LINEAR, GL_CLAMP_TO_EDGE, &textureGradient))
				throw std::runtime_error("error creating gradient texture");
		}

		GLuint textureOccupancyMap;
		if (!TextureHelper::CreateTexture3D(GL_R8UI, sizeOccupancyMap[0], sizeOccupancyMap[1], sizeOccupancyMap[2], GL_RED_INTEGER, GL_UNSIGNED_BYTE, NULL, GL_NEAREST, GL_CLAMP_TO_BORDER, &textureOccupancyMap))
			throw std::runtime_error("error creating occupany map texture");

		GLuint textureDistanceMap;
		if (!TextureHelper::CreateTexture3D(GL_R8UI, sizeOccupancyMap[0], sizeOccupancyMap[1], sizeOccupancyMap[2], GL_RED_INTEGER, GL_UNSIGNED_BYTE, NULL, GL_NEAREST, GL_CLAMP_TO_BORDER, &textureDistanceMap))
			throw std::runtime_error("error creating distance map texture");

		GLuint textureSwapArray[2];
		if (!TextureHelper::CreateTexture3D(GL_R8UI, sizeOccupancyMap[0], sizeOccupancyMap[1], sizeOccupancyMap[2], GL_RED_INTEGER, GL_UNSIGNED_BYTE, NULL, GL_NEAREST, GL_CLAMP_TO_BORDER, textureSwapArray, 2))
			throw std::runtime_error("error creating swap textures");

		GLuint textureSwap = textureSwapArray[0];

		unique_ptr<GLuint> texturePOMs(new GLuint[numIntensityPartitions]);
		if (!TextureHelper::CreateTexture3D(GL_R8UI, sizeOccupancyMap[0], sizeOccupancyMap[1], sizeOccupancyMap[2], GL_RED_INTEGER, GL_UNSIGNED_BYTE, NULL, GL_NEAREST, GL_CLAMP_TO_BORDER, texturePOMs.get(), numIntensityPartitions))
			throw std::runtime_error("error creating partitioned occupancy map textures");

		unique_ptr<GLuint> texturePDMs(new GLuint[numIntensityPartitions]);
		if (!TextureHelper::CreateTexture3D(GL_R8UI, sizeOccupancyMap[0], sizeOccupancyMap[1], sizeOccupancyMap[2], GL_RED_INTEGER, GL_UNSIGNED_BYTE, NULL, GL_NEAREST, GL_CLAMP_TO_BORDER, texturePDMs.get(), numIntensityPartitions))
			throw std::runtime_error("error creating partitioned distance map textures");

		if (!TextureHelper::CreateTexture2D(GL_R8, window_width, window_height, GL_RED, GL_UNSIGNED_BYTE, imgNoise.GetBufferAsVoid(), GL_NEAREST, GL_REPEAT, &textureNoise))
			throw std::runtime_error("error creating noise texture");

		GLuint textureITF;
		if (!TextureHelper::CreateTexture2D(GL_RGBA8, 256, 256, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, (void*)itfData, GL_NEAREST, GL_CLAMP_TO_EDGE, &textureITF))
			throw std::runtime_error("error creating ITF texture");


		// execute compute shaders
		if (precompute_gradient)
		{
			ComputeGradientMap(*computeshader_compute_gradient_map, textureVolume, voxelOutputDataType, textureGradient, volumeSize, 1.0);
			shaderParams.gradientType = 1;
		}

		ComputeOccupancyMap(*computeshader_compute_occupancy_map, textureVolume, voxelOutputDataType, textureGradient, textureITF, textureOccupancyMap, sizeOccupancyMap, reductionBlockSize, false);
		void* function = static_cast<void (*)(const ComputeShaderType&, const GLuint&, const VoxelOutputDataType, const GLuint&, const GLuint&, const GLuint&, const int*, const unsigned int&, const bool&)>(&ComputeOccupancyMap);
		DoFunctionTimings<void, const ComputeShaderType&, const GLuint&, const VoxelOutputDataType, const GLuint&, const GLuint&, const GLuint&, const int*, const unsigned int&, const bool&>(function, "occupancy map", numTimingIterations, *computeshader_compute_occupancy_map, textureVolume, voxelOutputDataType, textureGradient, textureITF, textureOccupancyMap, sizeOccupancyMap, reductionBlockSize, false);

		ComputeDistanceMap(computeshader_compute_distance_map, textureOccupancyMap, textureDistanceMap, textureSwap, sizeOccupancyMap);
		function = static_cast<void (*)(const ComputeShaderType&, const GLuint&, const GLuint&, const GLuint&, const int*)>(&ComputeDistanceMap);
		DoFunctionTimings<void, const ComputeShaderType&, const GLuint&, const GLuint&, const GLuint&, const int*>(function, "distance map", numTimingIterations, computeshader_compute_distance_map, textureOccupancyMap, textureDistanceMap, textureSwap, sizeOccupancyMap);

		std::vector<float> intensitiesPartitions_minValues(numIntensityPartitions);
		for (int i = 0; i < intensitiesPartitions_minValues.size(); i++)
			intensitiesPartitions_minValues[i] = ceil(i * (float)255 / (float)numIntensityPartitions);
		std::vector<float> intensitiesPartitions_maxValues(numIntensityPartitions);
		for (int i = 0; i < intensitiesPartitions_maxValues.size(); i++)
			intensitiesPartitions_maxValues[i] = floor((i + 1) * (float)255 / (float)numIntensityPartitions);

		if (useSpecialPartitionForZeroIntensity)
		{
			int numFreeIntensityPartitions = numIntensityPartitions - (int)(useSpecialPartitionForZeroIntensity); // reduce number of free intensity partitions in case of zero intensity special partition use
			for (int i = 1; i < intensitiesPartitions_minValues.size(); i++)
				intensitiesPartitions_minValues[i] = ceil((i - 1) * (float)255 / (float)numFreeIntensityPartitions);
			for (int i = 1; i < intensitiesPartitions_maxValues.size(); i++)
				intensitiesPartitions_maxValues[i] = floor((i) * (float)255 / (float)numFreeIntensityPartitions);
			intensitiesPartitions_minValues[0] = 0;
			intensitiesPartitions_maxValues[0] = 0;
			if (numIntensityPartitions > 1)
				intensitiesPartitions_minValues[1] = 1;
		}

		std::vector<bool> PDMsToUse(numIntensityPartitions);
		std::fill(PDMsToUse.begin(), PDMsToUse.end(), false);

		for (int i = 0; i < numIntensityPartitions; i++)
		{
			bool useMap = false;

			for (uint32_t intensity = intensitiesPartitions_minValues[i]; intensity <= intensitiesPartitions_maxValues[i]; intensity++)
			{
				useMap |= itfInterpolated.GetMappedAlphaValue(intensity) > 0;
			}
			PDMsToUse[i] = useMap;
		}

		ComputePOM(*computeshader_compute_POM, textureVolume, voxelOutputDataType, textureGradient, texturePOMs.get()[0], sizeOccupancyMap, reductionBlockSize, intensitiesPartitions_minValues[0], intensitiesPartitions_maxValues[0]);
		glFlush();
		glFinish();
		StartCounter();
		for (int i = 0; i < numIntensityPartitions; i++)
			ComputePOM(*computeshader_compute_POM, textureVolume, voxelOutputDataType, textureGradient, texturePOMs.get()[i], sizeOccupancyMap, reductionBlockSize, intensitiesPartitions_minValues[i], intensitiesPartitions_maxValues[i]);
		glFlush();
		glFinish();
		std::cout << "partitioned occupancy maps computation took " << GetCounter() << " msecs\n";

		ComputeDistanceMap(computeshader_compute_distance_map, texturePOMs.get()[0], texturePDMs.get()[0], textureSwap, sizeOccupancyMap);
		glFlush();
		glFinish();
		StartCounter();
		for (int i = 0; i < numIntensityPartitions; i++)
			ComputeDistanceMap(computeshader_compute_distance_map, texturePOMs.get()[i], texturePDMs.get()[i], textureSwap, sizeOccupancyMap);
		glFlush();
		glFinish();
		std::cout << "partitioned distance maps computation took " << GetCounter() << " msecs\n";

		if (!disablePartitionedDistanceMapOptimization)
		{
			CombinePDMs(computeshader_combine_PDMs, texturePDMs.get(), textureSwapArray, textureDistanceMap, sizeOccupancyMap, numIntensityPartitions, PDMsToUse, false);
			function = static_cast<void (*)(const ComputeShaderType&, const GLuint*, const GLuint*, const GLuint&, int*, const unsigned int&, const std::vector<bool>, const bool&)>(&CombinePDMs);
			DoFunctionTimings<void, const ComputeShaderType&, const GLuint*, const GLuint*, const GLuint&, int*, const unsigned int&, const std::vector<bool>, const bool&>(function, "combining distance maps", numTimingIterations, computeshader_combine_PDMs, texturePDMs.get(), textureSwapArray, textureDistanceMap, sizeOccupancyMap, numIntensityPartitions, PDMsToUse, false);
		}

		projection = glm::perspective((float)glm::radians(fov), (float)window_width / (float)window_height, 0.3f, 100.0f);

		float physical_size[3] = { (float)volumeSize[0] * (float)voxelSize[0], (float)volumeSize[1] * (float)voxelSize[1], (float)volumeSize[2] * (float)voxelSize[2] };
		glm::mat4 scaling_cube = glm::scale(glm::mat4(1.0f), glm::vec3(physical_size[0], physical_size[1], physical_size[2]));

		glm::mat4 translation_cube = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));

		glm::mat4 volumeRotationFromSitkDirection = glm::mat4(direction[0], direction[1], direction[2], 0.0, direction[3], direction[4], direction[5], 0.0, direction[6], direction[7], direction[8], 0.0, 0.0, 0.0, 0.0, 1.0);

		glm::mat4 initial_rotation_x = glm::rotate(glm::mat4(1.0f), initialVolumeRotationX / 180.0f * PI, glm::vec3(1.0f, 0.0f, 0.0f));
		glm::mat4 initial_rotation_y = glm::rotate(glm::mat4(1.0f), initialVolumeRotationY / 180.0f * PI, glm::vec3(0.0f, 1.0f, 0.0f));
		glm::mat4 initial_rotation_z = glm::rotate(glm::mat4(1.0f), initialVolumeRotationZ / 180.0f * PI, glm::vec3(0.0f, 0.0f, 1.0f));

		glm::mat4 transformation_initial = translation_cube * initial_rotation_z * initial_rotation_y * initial_rotation_x * volumeRotationFromSitkDirection * scaling_cube;

		glm::mat4 custom_rotation_y = glm::mat4(1.0f);

		if (useDemoAndTimingsMode)
		{
			cam.SetCameraPos(glm::vec3(0.0f, demoAndTimingsModeCameraHeight, demoAndTimingsModeCameraDistance));
			cam.SetCameraPitch(demoAndTimingsModeCameraPitch);
		}
		glfwSetWindowSizeCallback(window, WindowResizeCallback);

		float samplingFactor = 1.0f;
		int itf_gamma_correction = (int)itfInterpolated.getApplyGammaCorrection();
		float gamma = itfInterpolated.getGamma();
		glm::vec3 lightDirection = glm::vec3(-0.5f, -1.0f, 0.0f);

		unsigned int ubo_vert_block;
		glGenBuffers(1, &ubo_vert_block);
		glBindBuffer(GL_UNIFORM_BUFFER, ubo_vert_block);
		glBufferData(GL_UNIFORM_BUFFER, 64, NULL, GL_STATIC_DRAW);
		glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_vert_block);
		glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_vert_block, 0, 64);

		unsigned int ubo_frag_block;
		glGenBuffers(1, &ubo_frag_block);
		glBindBuffer(GL_UNIFORM_BUFFER, ubo_frag_block);
		glBufferData(GL_UNIFORM_BUFFER, 268, NULL, GL_STATIC_DRAW);
		glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo_frag_block);
		glBindBufferRange(GL_UNIFORM_BUFFER, 1, ubo_frag_block, 0, 268);

		glBindBuffer(GL_UNIFORM_BUFFER, 0);

		unsigned int loc_shader0_vertexPosition = 0;

		loc_shader0_vertexPosition = glGetAttribLocation(shader_render_volume.ID, "vertexPosition");
		glEnableVertexAttribArray(loc_shader0_vertexPosition);
		glVertexAttribPointer(loc_shader0_vertexPosition, 3, GL_FLOAT, GL_FALSE, 12, 0);
		shader_render_volume.use();

		double currentTime = glfwGetTime();

		double lastTime = glfwGetTime();
		unsigned int numFramesInCurrentSecond = 0;

		int frameCountSumRenderingTiming = 0;
		double lastFrameTime = 0.0;
		double firstMeasuredFrameTimeStart = 0.0;

		while (!glfwWindowShouldClose(window))
		{
			glfwPollEvents();

			DoMovement();

			glm::mat4 view = cam.GetViewMatrix();

			float* data = glm::value_ptr(view);

			if (useDemoAndTimingsMode)
				custom_rotation_y = glm::rotate(glm::mat4(1.0f), (float)glfwGetTime() * (rotationsPerTimingPeriod / (float)timingPeriodInSec) * 2.0f * PI, glm::vec3(0.0f, 1.0f, 0.0f));

			glm::mat4 transformation_cube = translation_cube * custom_rotation_y * transformation_initial;
			glm::mat4 modelViewProjection = projection * view * transformation_cube;

			glBindVertexArray(renderGeometry.VAO_cube);

			// render back faces of cube and do volume rendering
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			glCullFace(GL_FRONT);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_3D, textureVolume);

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_3D, textureGradient);

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, textureNoise);

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, textureITF);

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_3D, textureOccupancyMap);

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_3D, textureDistanceMap);

			glm::mat4x4 inverseModelView = glm::inverse(view * transformation_cube);

			glBindBuffer(GL_UNIFORM_BUFFER, ubo_vert_block);
			glBufferSubData(GL_UNIFORM_BUFFER, 0, 64, glm::value_ptr(modelViewProjection));

			glBindBuffer(GL_UNIFORM_BUFFER, ubo_frag_block);
			glBufferSubData(GL_UNIFORM_BUFFER, 0, 64, glm::value_ptr(transformation_cube));
			glBufferSubData(GL_UNIFORM_BUFFER, 64, 64, glm::value_ptr(view));
			glBufferSubData(GL_UNIFORM_BUFFER, 128, 64, glm::value_ptr(inverseModelView));
			glBufferSubData(GL_UNIFORM_BUFFER, 192, 4, &samplingFactor);
			glBufferSubData(GL_UNIFORM_BUFFER, 196, 4, &shaderParams.renderContentType);
			glBufferSubData(GL_UNIFORM_BUFFER, 200, 4, &shaderParams.ditheringMode);
			glBufferSubData(GL_UNIFORM_BUFFER, 204, 4, &shaderParams.earlyRayTermination);
			glBufferSubData(GL_UNIFORM_BUFFER, 208, 4, &shaderParams.essType);
			glBufferSubData(GL_UNIFORM_BUFFER, 212, 4, &shaderParams.gradientType);
			glBufferSubData(GL_UNIFORM_BUFFER, 216, 4, &itf_gamma_correction);
			glBufferSubData(GL_UNIFORM_BUFFER, 220, 4, &gamma);
			glBufferSubData(GL_UNIFORM_BUFFER, 224, 4, &enableTricubicFiltering);
			glBufferSubData(GL_UNIFORM_BUFFER, 228, 4, &shaderParams.reuseNormals);
			glBufferSubData(GL_UNIFORM_BUFFER, 232, 4, &shaderParams.lightshadingFactor);
			glBufferSubData(GL_UNIFORM_BUFFER, 240, 4, &shaderParams.lightDirComponentX);
			glBufferSubData(GL_UNIFORM_BUFFER, 244, 4, &shaderParams.lightDirComponentY);
			glBufferSubData(GL_UNIFORM_BUFFER, 248, 4, &shaderParams.lightDirComponentZ);
			glBufferSubData(GL_UNIFORM_BUFFER, 252, 4, &shaderParams.lightAmbientFactor);
			glBufferSubData(GL_UNIFORM_BUFFER, 256, 4, &shaderParams.lightDiffuseFactor);
			glBufferSubData(GL_UNIFORM_BUFFER, 260, 4, &shaderParams.lightSpecularFactor);
			glBufferSubData(GL_UNIFORM_BUFFER, 264, 4, &shaderParams.lightSpecularShininess);

			glBindBuffer(GL_UNIFORM_BUFFER, 0);

			glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);

			glfwSwapBuffers(window);

			currentTime = glfwGetTime();

			if (useDemoAndTimingsMode)
			{
				if (lastFrameTime != 0.0)
				{
					frameCountSumRenderingTiming++;
					if (currentTime - firstMeasuredFrameTimeStart >= timingPeriodInSec)
					{
						std::cout << "mean rendering time for a single frame over " << timingPeriodInSec << " seconds was " << ((timingPeriodInSec * 1000.0) / (double)frameCountSumRenderingTiming) << " msecs\n";
						frameCountSumRenderingTiming = 0.0;
						firstMeasuredFrameTimeStart = currentTime;
					}
				}
				else
				{
					firstMeasuredFrameTimeStart = currentTime;
				}
			}

			if (itfInterpolate > 0)
			{
				if (itfInterpolate == 1)
				{
					t_itfInterpolation -= ((currentTime - lastFrameTime)) / itfBlendTime;
				}
				if (itfInterpolate == 2)
				{
					t_itfInterpolation += ((currentTime - lastFrameTime)) / itfBlendTime;
				}
				t_itfInterpolation = std::clamp(t_itfInterpolation, 0.0f, 1.0f);
				itfInterpolated.SetFromInterpolatedITFs(itf, itf2, t_itfInterpolation);
				itfData = itfInterpolated.GetData();
				itf_gamma_correction = (int)itfInterpolated.getApplyGammaCorrection();
				gamma = itfInterpolated.getGamma();

				GLfloat borderColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
				glBindTexture(GL_TEXTURE_2D, textureITF);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 256, 256, 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, itfData);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
				glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

				if (!disablePartitionedDistanceMapOptimization)
				{
					std::vector<bool> PDMsToUse(numIntensityPartitions);
					std::fill(PDMsToUse.begin(), PDMsToUse.end(), false);
					for (unsigned int i = 0; i < numIntensityPartitions; i++)
					{
						bool useMap = false;

						for (uint32_t intensity = intensitiesPartitions_minValues[i]; intensity <= intensitiesPartitions_maxValues[i]; intensity++)
						{
							useMap |= itfInterpolated.GetMappedAlphaValue(intensity) > 0;
						}
						PDMsToUse[i] = useMap;
					}
					CombinePDMs(computeshader_combine_PDMs, texturePDMs.get(), textureSwapArray, textureDistanceMap, sizeOccupancyMap, numIntensityPartitions, PDMsToUse, false);
				}
				else
				{
					ComputeDistanceMap(computeshader_compute_distance_map, textureOccupancyMap, textureDistanceMap, textureSwap, sizeOccupancyMap);
				}
				// reset ubo for binding 0, since compute shader used binding for its own ubo
				glBindBuffer(GL_UNIFORM_BUFFER, ubo_vert_block);
				glBufferData(GL_UNIFORM_BUFFER, 64, NULL, GL_STATIC_DRAW);
				glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_vert_block);
				glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_vert_block, 0, 64);
				glBindBuffer(GL_UNIFORM_BUFFER, 0);

				shader_render_volume.use();

				if (t_itfInterpolation >= 1.0f)
				{
					itfInterpolate = 0;
					t_itfInterpolation = 1.0f;
				}
				if (t_itfInterpolation <= 0.0f)
				{
					itfInterpolate = 0;
					t_itfInterpolation = 0.0f;
				}
			}

			lastFrameTime = currentTime;

			numFramesInCurrentSecond += 1;
			if (currentTime - lastTime >= 1.0)
			{
				string windowTitle = string("frame drawing took ") + std::format("{:.3f}", 1000.0 / (double)numFramesInCurrentSecond) + string(" ms") +
					string(" rendering mode: ") + std::format("{}", shaderParams.renderContentType) +
					string(" dithering mode: ") + std::format("{}", shaderParams.ditheringMode) +
					string(" early ray termination: ") + std::format("{}", shaderParams.earlyRayTermination) +
					string(" ess type: ") + std::format("{}", shaderParams.essType) +
					string(" light shading: ") + std::format("{}", shaderParams.lightshadingFactor) +
					string(" precompute gradient: ") + std::format("{}", shaderParams.gradientType) +
					string(" tricubic filtering: ") + std::format("{}", enableTricubicFiltering) +
					string(" reuse normals: ") + std::format("{}", shaderParams.reuseNormals);

				glfwSetWindowTitle(window, windowTitle.c_str());

				numFramesInCurrentSecond = 0;
				lastTime = glfwGetTime();
			}
		}

		TextureHelper::DeleteTexture(textureVolume);
		TextureHelper::DeleteTexture(textureGradient);
		TextureHelper::DeleteTexture(textureDistanceMap);
		TextureHelper::DeleteTextures(2, textureSwapArray);
		TextureHelper::DeleteTextures(numIntensityPartitions, texturePOMs.get());
		TextureHelper::DeleteTextures(numIntensityPartitions, texturePDMs.get());

		glfwDestroyWindow(window);

		glfwTerminate();
	}
	catch (const std::exception& err)
	{
		std::cerr << "error occured:\n" << err.what() << "\n";
		return 1;
	}

	return 0;
}
