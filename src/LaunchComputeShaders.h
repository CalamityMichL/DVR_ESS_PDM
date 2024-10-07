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

#pragma once

#include "./common/ComputeShaderSpirv.h"

#include <glad/gl.h>
#include <math.h>
#include <vector>

#include <glm/vec4.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef enum VoxelDataType
{
    DTYPE_UNORM = 0,
    DTYPE_INT16 = 1,
    DTYPE_UINT16 = 2,
    DTYPE_UINT8 = 3
} VoxelDataType;

typedef enum VoxelOutputDataType
{
    DTYPE_UNORM8 = 0,
    DTYPE_UNORM16 = 1
} VoxelOutputDataType;

void _ComputeConvertVolumeToSNORM(const GLuint& textureVolume, const VoxelDataType datatype, const VoxelOutputDataType outputDatatype, const VoxelOutputDataType outputDataType, const float minIntensityValue, const float maxIntensityValue, const GLuint& textureOutputVolume, std::vector<unsigned int>& volumeSize);
void ComputeConvertVolumeToSNORM(const ComputeShaderSpirv& shader_compute_convert_volume_to_SNORM, const GLuint& textureVolume, const VoxelDataType datatype, const VoxelOutputDataType outputDatatype, const float minIntensityValue, const float maxIntensityValue, const GLuint& textureOutputVolume, std::vector<unsigned int>& volumeSize);

void _ComputeGradientMap(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, std::vector<unsigned int>& volumeSize, const float& grad_magnitude_modifier);
void ComputeGradientMap(const ComputeShaderSpirv& shader_compute_gradient, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, std::vector<unsigned int>& volumeSize, const float& grad_magnitude_modifier);

void _ComputeOccupancyMap(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& textureITF, const GLuint& textureOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const bool& useITF2D);
void ComputeOccupancyMap(const ComputeShaderSpirv& shader_compute_occupancy, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& textureITF, const GLuint& textureOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const bool& useITF2D);

void _ComputePOM(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& texturePartitionedOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const float& partition_intensity_min, const float& partition_intensity_max);
void ComputePOM(const ComputeShaderSpirv& shader_compute_partitioned_occupancy, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& texturePartitionedOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const float& partition_intensity_min, const float& partition_intensity_max);

void _ComputeDistanceMap(const GLuint& textureOccupancyMap, const GLuint& textureDistanceMap, const GLuint& textureSwap, const int* sizeOccupancyMap);
void ComputeDistanceMap(const ComputeShaderSpirv& shader_compute_distance, const GLuint& textureOccupancyMap, const GLuint& textureDistanceMap, const GLuint& textureSwap, const int* sizeOccupancyMap);

void _CombinePDMs(const GLuint* texturePartitionedDistanceMaps, const GLuint* textureSwapArray, const GLuint& textureIntermediateDistanceMapOutput, int* sizeOccupancyMap, const unsigned int& numIntensityPartitions, const std::vector<bool> partitionedDistanceMapsToUse, const bool& debugOut);
void CombinePDMs(const ComputeShaderSpirv& shader_combine_distance_maps, const GLuint* texturePartitionedDistanceMaps, const GLuint* textureSwapArray, const GLuint& textureIntermediateDistanceMapOutput, int* sizeOccupancyMap, const unsigned int& numIntensityPartitions, const std::vector<bool> partitionedDistanceMapsToUse, const bool& debugOut);

void _ComputeVoxelHistogram(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureHistogram, std::vector<unsigned int>& volumeSize);
void ComputeVoxelHistogram(const ComputeShaderSpirv& shader_compute_histogram, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureHistogram, std::vector<unsigned int>& volumeSize);
