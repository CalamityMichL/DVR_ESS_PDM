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
 *
 * code derived from: https://github.com/JoeyDeVries/LearnOpenGL/blob/master/includes/learnopengl/shader.h
 *
 */

#pragma once

#include <glad/gl.h>
#include <glm/glm.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class ComputeShaderSpirv
{
public:
    unsigned int ID;
    // constructor generates the shader on the fly
    // ------------------------------------------------------------------------
    ComputeShaderSpirv(const char* computeBinaryPath)
    {
        // 1. retrieve the compute source code from computeBinaryPath
        std::unique_ptr<uint8_t[]> bufferComputeShader;
        size_t sizeComputeBuffer;       
        std::ifstream computeShaderFile;        
        // ensure ifstream objects can throw exceptions:
        computeShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try
        {
            computeShaderFile = std::ifstream(computeBinaryPath, std::ios::in | std::ios::binary | std::ios::ate);
            if (!computeShaderFile)
                throw std::ifstream::failure("input spirv compute shader file not found (" + std::string(computeBinaryPath) + ")");
            sizeComputeBuffer = computeShaderFile.tellg();
            computeShaderFile.seekg(0, std::ios::beg);
            bufferComputeShader = std::make_unique<uint8_t[]>(sizeComputeBuffer);
            computeShaderFile.read((char*)bufferComputeShader.get(), sizeComputeBuffer);
            computeShaderFile.close();
        }
        catch (std::ifstream::failure& e)
        {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << e.what() << std::endl;
            throw(e);
        }

        // 2. compile compute shader
        unsigned int compute;        
        compute = glCreateShader(GL_COMPUTE_SHADER);
        glShaderBinary(1, &compute, GL_SHADER_BINARY_FORMAT_SPIR_V_ARB, bufferComputeShader.get(), (GLsizei)sizeComputeBuffer);
        glSpecializeShader(compute, "main", 0, 0, 0);
        CheckCompileErrors(compute, "COMPUTE");

        // shader Program
        ID = glCreateProgram();
        glAttachShader(ID, compute);
        glLinkProgram(ID);
        CheckCompileErrors(ID, "PROGRAM");
        // delete the shaders as they're linked into our program now and no longer necessary
        glDeleteShader(compute);

    }

    ~ComputeShaderSpirv() {
        glDeleteProgram(ID);
    }

    // activate the shader
    // ------------------------------------------------------------------------
    void use() const
    {
        glUseProgram(ID);
    }

private:
    // utility function for checking shader compilation/linking errors.
    // ------------------------------------------------------------------------
    void CheckCompileErrors(GLuint shader, std::string type)
    {
        GLint success;
        GLchar infoLog[1024];
        if (type != "PROGRAM")
        {
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success)
            {
                glGetShaderInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
        else
        {
            glGetProgramiv(shader, GL_LINK_STATUS, &success);
            if (!success)
            {
                glGetProgramInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
    }
};
