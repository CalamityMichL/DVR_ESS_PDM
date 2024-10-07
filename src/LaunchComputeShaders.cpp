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

#include "LaunchComputeShaders.h"

#include "./common/ComputeShaderSpirv.h"

#include <glad/gl.h>
#include <math.h>
#include <vector>

#include <glm/vec4.hpp>
#include <glm/gtc/type_ptr.hpp>

void _ComputeConvertVolumeToSNORM(const GLuint& textureVolume, const VoxelDataType datatype, const VoxelOutputDataType outputDatatype, const float minIntensityValue, const float maxIntensityValue, const GLuint& textureOutputVolume, std::vector<unsigned int>& volumeSize)
{
    if (datatype == VoxelDataType::DTYPE_INT16)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R16I);
    else if (datatype == VoxelDataType::DTYPE_UINT16)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R16UI);
    else if (datatype == VoxelDataType::DTYPE_UINT8)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8UI);

    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM8)
        glBindImageTexture(1, textureOutputVolume, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8);
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM16)
        glBindImageTexture(1, textureOutputVolume, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R16);

    unsigned int ubo_block;
    glGenBuffers(1, &ubo_block);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_block);
    glBufferData(GL_UNIFORM_BUFFER, 8, NULL, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_block);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_block, 0, 8);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 4, &minIntensityValue);
    glBufferSubData(GL_UNIFORM_BUFFER, 4, 4, &maxIntensityValue);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)volumeSize[0] / 8), (unsigned int)ceil((float)volumeSize[1] / 8), (unsigned int)ceil((float)volumeSize[2] / 8) };
    glDispatchCompute(numWorkGroups[0], numWorkGroups[1], numWorkGroups[2]);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glDeleteBuffers(1, &ubo_block);
}

void ComputeConvertVolumeToSNORM(const ComputeShaderSpirv& shader_compute_convert_volume_to_SNORM, const GLuint& textureVolume, const VoxelDataType datatype, const VoxelOutputDataType outputDatatype, const float minIntensityValue, const float maxIntensityValue, const GLuint& textureOutputVolume, std::vector<unsigned int>& volumeSize)
{
    shader_compute_convert_volume_to_SNORM.use();
    _ComputeConvertVolumeToSNORM(textureVolume, datatype, outputDatatype, minIntensityValue, maxIntensityValue, textureOutputVolume, volumeSize);
}

void _ComputeGradientMap(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, std::vector<unsigned int>& volumeSize, const float& grad_magnitude_modifier)
{
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM8)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8);
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM16)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R16);

    glBindImageTexture(1, textureGradient, 0, GL_TRUE, 0, GL_READ_WRITE, GL_RGBA8_SNORM);

    unsigned int ubo_block;
    glGenBuffers(1, &ubo_block);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_block);
    glBufferData(GL_UNIFORM_BUFFER, 4, NULL, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_block);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_block, 0, 4);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 4, &grad_magnitude_modifier);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)volumeSize[0] / 8), (unsigned int)ceil((float)volumeSize[1] / 8), (unsigned int)ceil((float)volumeSize[2] / 8) };
    glDispatchCompute(numWorkGroups[0], numWorkGroups[1], numWorkGroups[2]);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glDeleteBuffers(1, &ubo_block);
}
void ComputeGradientMap(const ComputeShaderSpirv& shader_compute_gradient, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, std::vector<unsigned int>& volumeSize, const float& grad_magnitude_modifier)
{
    shader_compute_gradient.use();
    _ComputeGradientMap(textureVolume, outputDatatype, textureGradient, volumeSize, grad_magnitude_modifier);
}

void _ComputeOccupancyMap(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& textureITF, const GLuint& textureOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const bool& useITF2D)
{
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM8)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8);
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM16)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R16);

    glBindImageTexture(1, textureGradient, 0, GL_TRUE, 0, GL_READ_ONLY, GL_RGBA8);
    glBindImageTexture(2, textureITF, 0, GL_TRUE, 0, GL_READ_ONLY, GL_RGBA8);
    glBindImageTexture(3, textureOccupancyMap, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);

    unsigned int ubo_block;
    glGenBuffers(1, &ubo_block);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_block);
    glBufferData(GL_UNIFORM_BUFFER, 20, NULL, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_block);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_block, 0, 20);
    glm::vec4 vecReductionBlockSize = glm::vec4((float)reductionBlockSize, (float)reductionBlockSize, (float)reductionBlockSize, 0.0f);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 16, glm::value_ptr(vecReductionBlockSize));
    int useITF2D_int = (int)useITF2D;
    glBufferSubData(GL_UNIFORM_BUFFER, 16, 4, &useITF2D_int);

    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)sizeOccupancyMap[0] / 8), (unsigned int)ceil((float)sizeOccupancyMap[1] / 8), (unsigned int)ceil((float)sizeOccupancyMap[2] / 8) };
    glDispatchCompute(numWorkGroups[0], numWorkGroups[1], numWorkGroups[2]);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glDeleteBuffers(1, &ubo_block);
}
void ComputeOccupancyMap(const ComputeShaderSpirv& shader_compute_occupancy, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& textureITF, const GLuint& textureOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const bool& useITF2D)
{
    shader_compute_occupancy.use();
    _ComputeOccupancyMap(textureVolume, outputDatatype, textureGradient, textureITF, textureOccupancyMap, sizeOccupancyMap, reductionBlockSize, useITF2D);
}


void _ComputePOM(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& texturePartitionedOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const float& partition_intensity_min = 0, const float& partition_intensity_max = 255)
{
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM8)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8);
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM16)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R16);

    glBindImageTexture(1, textureGradient, 0, GL_TRUE, 0, GL_READ_ONLY, GL_RGBA8);
    // glBindImageTexture(2, textureITF, 0, GL_TRUE, 0, GL_READ_ONLY, GL_RGBA8);
    glBindImageTexture(2, texturePartitionedOccupancyMap, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);

    unsigned int ubo_block;
    glGenBuffers(1, &ubo_block);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_block);
    glBufferData(GL_UNIFORM_BUFFER, 24, NULL, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_block);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_block, 0, 24);
    glm::vec4 vecReductionBlockSize = glm::vec4((float)reductionBlockSize, (float)reductionBlockSize, (float)reductionBlockSize, 0.0f);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 16, glm::value_ptr(vecReductionBlockSize));
    glBufferSubData(GL_UNIFORM_BUFFER, 16, 4, &partition_intensity_min);
    glBufferSubData(GL_UNIFORM_BUFFER, 20, 4, &partition_intensity_max);

    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)sizeOccupancyMap[0] / 8), (unsigned int)ceil((float)sizeOccupancyMap[1] / 8), (unsigned int)ceil((float)sizeOccupancyMap[2] / 8) };
    glDispatchCompute(numWorkGroups[0], numWorkGroups[1], numWorkGroups[2]);
    glMemoryBarrier(GL_ALL_BARRIER_BITS); // (GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    glDeleteBuffers(1, &ubo_block);
}
void ComputePOM(const ComputeShaderSpirv& shader_compute_partitioned_occupancy, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureGradient, const GLuint& texturePartitionedOccupancyMap, const int* sizeOccupancyMap, const unsigned int& reductionBlockSize, const float& partition_intensity_min = 0, const float& partition_intensity_max = 255)
{
    shader_compute_partitioned_occupancy.use();
    _ComputePOM(textureVolume, outputDatatype, textureGradient, texturePartitionedOccupancyMap, sizeOccupancyMap, reductionBlockSize, partition_intensity_min, partition_intensity_max);
}

void _ComputeDistanceMap(const GLuint& textureOccupancyMap, const GLuint& textureDistanceMap, const GLuint& textureSwap, const int* sizeOccupancyMap)
{
    glBindImageTexture(0, textureDistanceMap, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);
    glBindImageTexture(1, textureOccupancyMap, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);

    unsigned int ubo_block;
    glGenBuffers(1, &ubo_block);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_block);
    glBufferData(GL_UNIFORM_BUFFER, 4, NULL, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_block);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_block, 0, 4);


    unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)sizeOccupancyMap[0] / 8), (unsigned int)ceil((float)sizeOccupancyMap[1] / 8), (unsigned int)ceil((float)sizeOccupancyMap[2] / 8) };

    int stage = 0;
    // Dispatch 1st stage
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 4, &stage);
    glDispatchCompute(numWorkGroups[1], numWorkGroups[2], 1);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    // Dispatch 2nd stage
    glBindImageTexture(1, textureSwap, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);
    stage = 1;
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 4, &stage);
    glDispatchCompute(numWorkGroups[0], numWorkGroups[2], 1);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    // Dispatch 3rd stage
    stage = 2;
    glBufferSubData(GL_UNIFORM_BUFFER, 0, 4, &stage);
    glDispatchCompute(numWorkGroups[0], numWorkGroups[1], 1);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    glDeleteBuffers(1, &ubo_block);
}
void ComputeDistanceMap(const ComputeShaderSpirv& shader_compute_distance, const GLuint& textureOccupancyMap, const GLuint& textureDistanceMap, const GLuint& textureSwap, const int* sizeOccupancyMap)
{
    shader_compute_distance.use();
    _ComputeDistanceMap(textureOccupancyMap, textureDistanceMap, textureSwap, sizeOccupancyMap);
}

void _CombinePDMs(const GLuint* texturePartitionedDistanceMaps, const GLuint* textureSwapArray, const GLuint& textureIntermediateDistanceMapOutput, int* sizeOccupancyMap, const unsigned int& numIntensityPartitions, const std::vector<bool> partitionedDistanceMapsToUse, const bool& debugOut = false)
{
    using namespace std;

    glBindImageTexture(0, textureSwapArray[0], 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);
    glBindImageTexture(1, textureSwapArray[1], 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);
    unsigned int numMapsToCombinePerLaunch = 6;
    unsigned int numMapsToCombineTotal = 0;
    for (int i = 0; i < partitionedDistanceMapsToUse.size(); i++)
        numMapsToCombineTotal += partitionedDistanceMapsToUse[i] == true;
    int numIterations = 0;
    if (numMapsToCombineTotal > 0)
        numIterations = (int)ceil((float)(numMapsToCombineTotal) / (float)(numMapsToCombinePerLaunch));
    unsigned int idxMapProcessed = 0;
    unsigned int countSkippedMaps = 0;
    int j = 0;
    for (; j < numIterations; j++)
    {
        glBindImageTexture(0, textureSwapArray[j % 2], 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);
        glBindImageTexture(1, textureSwapArray[(j + 1) % 2], 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8UI);
        unsigned int numMapsToUse = std::min(numMapsToCombinePerLaunch, numMapsToCombineTotal - j * numMapsToCombinePerLaunch);
        for (unsigned int i = 0; i < numMapsToUse; i++)
        {
            while (partitionedDistanceMapsToUse[idxMapProcessed] == false)
            {
                idxMapProcessed += 1;
                countSkippedMaps += 1;
            }
            glBindImageTexture(2 + i, texturePartitionedDistanceMaps[idxMapProcessed], 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8UI);
            idxMapProcessed += 1;
        }

        // glBindImageTexture(2, textureSubdividedDistanceMaps[1], 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8UI)

        unsigned int ubo_block;
        glGenBuffers(1, &ubo_block);
        glBindBuffer(GL_UNIFORM_BUFFER, ubo_block);
        glBufferData(GL_UNIFORM_BUFFER, 8, NULL, GL_STATIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_block);
        glBindBufferRange(GL_UNIFORM_BUFFER, 0, ubo_block, 0, 8);
        int useDistanceValueFromPreviousPass = int(j > 0);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, 4, &useDistanceValueFromPreviousPass);
        glBufferSubData(GL_UNIFORM_BUFFER, 4, 4, &numMapsToUse);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)sizeOccupancyMap[0] / 8), (unsigned int)ceil((float)sizeOccupancyMap[1] / 8), (unsigned int)ceil((float)sizeOccupancyMap[2] / 8) };
        glDispatchCompute(numWorkGroups[0], numWorkGroups[1], numWorkGroups[2]);
        glMemoryBarrier(GL_ALL_BARRIER_BITS);

        glDeleteBuffers(1, &ubo_block);
    }

    if (debugOut == true)
    {
        cout << "number of maps to be skipped according to skip struct: " << partitionedDistanceMapsToUse.size() - numMapsToCombineTotal << endl;
        cout << "number of skipped distance maps due to itf: " << countSkippedMaps + partitionedDistanceMapsToUse.size() - idxMapProcessed << endl;
    }

    glCopyImageSubData(textureSwapArray[(j) % 2], GL_TEXTURE_3D, 0, 0, 0, 0, textureIntermediateDistanceMapOutput, GL_TEXTURE_3D, 0, 0, 0, 0, int(sizeOccupancyMap[0]), int(sizeOccupancyMap[1]), int(sizeOccupancyMap[2]));
}
void CombinePDMs(const ComputeShaderSpirv& shader_combine_distance_maps, const GLuint* texturePartitionedDistanceMaps, const GLuint* textureSwapArray, const GLuint& textureIntermediateDistanceMapOutput, int* sizeOccupancyMap, const unsigned int& numIntensityPartitions, const std::vector<bool> partitionedDistanceMapsToUse, const bool& debugOut = false)
{
    shader_combine_distance_maps.use();
    _CombinePDMs(texturePartitionedDistanceMaps, textureSwapArray, textureIntermediateDistanceMapOutput, sizeOccupancyMap, numIntensityPartitions, partitionedDistanceMapsToUse, debugOut);
}


void _ComputeVoxelHistogram(const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureHistogram, std::vector<unsigned int>& volumeSize)
{
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM8)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R8);
    if (outputDatatype == VoxelOutputDataType::DTYPE_UNORM16)
        glBindImageTexture(0, textureVolume, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R16);
    glBindImageTexture(1, textureHistogram, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_R32UI);

    unsigned int numWorkGroups[3] = { (unsigned int)ceil((float)volumeSize[0] / 8), (unsigned int)ceil((float)volumeSize[1] / 8), (unsigned int)ceil((float)volumeSize[2] / 8) };
    glDispatchCompute(numWorkGroups[0], numWorkGroups[1], numWorkGroups[2]);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
}
void ComputeVoxelHistogram(const ComputeShaderSpirv& shader_compute_histogram, const GLuint& textureVolume, const VoxelOutputDataType outputDatatype, const GLuint& textureHistogram, std::vector<unsigned int>& volumeSize)
{
    shader_compute_histogram.use();
    _ComputeVoxelHistogram(textureVolume, outputDatatype, textureHistogram, volumeSize);
}
