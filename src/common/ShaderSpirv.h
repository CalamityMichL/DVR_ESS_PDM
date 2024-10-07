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

class ShaderSpirv
{
public:
    unsigned int ID;
    // constructor generates the shader on the fly
    // ------------------------------------------------------------------------
    ShaderSpirv(const char* vertexBinaryPath, const char* fragmentBinaryPath)
    {
        // 1. retrieve the vertex/fragment binary buffer from vertexBinaryPath and fragmentBinaryPath
        std::unique_ptr<uint8_t[]> bufferVertexShader;
        std::unique_ptr<uint8_t[]> bufferFragmentShader;
        size_t sizeVertexBuffer, sizeFragmentBuffer;
        std::ifstream vShaderFile;
        std::ifstream fShaderFile;
        // ensure ifstream objects can throw exceptions:
        vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try
        {
            // open files
            vShaderFile = std::ifstream(vertexBinaryPath, std::ios::in | std::ios::binary | std::ios::ate);
            if (!vShaderFile)
                throw std::ifstream::failure("input spirv vertex shader file not found (" + std::string(vertexBinaryPath) + ")");
            sizeVertexBuffer = vShaderFile.tellg();
            vShaderFile.seekg(0, std::ios::beg);
            bufferVertexShader = std::make_unique<uint8_t[]>(sizeVertexBuffer);
            vShaderFile.read((char*)bufferVertexShader.get(), sizeVertexBuffer);
            vShaderFile.close();

            fShaderFile = std::ifstream(fragmentBinaryPath, std::ios::in | std::ios::binary | std::ios::ate);
            if (!fShaderFile)
                throw std::ifstream::failure("input spirv fragment shader file not found (" + std::string(fragmentBinaryPath) + ")");
            sizeFragmentBuffer = fShaderFile.tellg();
            fShaderFile.seekg(0, std::ios::beg);
            bufferFragmentShader = std::make_unique<uint8_t[]>(sizeFragmentBuffer);
            fShaderFile.read((char*)bufferFragmentShader.get(), sizeFragmentBuffer);
            fShaderFile.close();
        }
        catch (std::ifstream::failure& e)
        {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << e.what() << std::endl;
            throw e;
        }

        // 2. compile shaders
        unsigned int vertex, fragment;
        // vertex shader
        vertex = glCreateShader(GL_VERTEX_SHADER);
        glShaderBinary(1, &vertex, GL_SHADER_BINARY_FORMAT_SPIR_V_ARB, bufferVertexShader.get(), (GLsizei)sizeVertexBuffer);
        glSpecializeShader(vertex, "main", 0, 0, 0);

        CheckCompileErrors(vertex, "VERTEX");
        // fragment Shader
        fragment = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderBinary(1, &fragment, GL_SHADER_BINARY_FORMAT_SPIR_V_ARB, bufferFragmentShader.get(), (GLsizei)sizeFragmentBuffer);
        glSpecializeShader(fragment, "main", 0, 0, 0);
        CheckCompileErrors(fragment, "FRAGMENT");
        // shader Program
        ID = glCreateProgram();
        glAttachShader(ID, vertex);
        glAttachShader(ID, fragment);
        glLinkProgram(ID);
        CheckCompileErrors(ID, "PROGRAM");
        // delete the shaders as they're linked into our program now and no longer necessary
        glDeleteShader(vertex);
        glDeleteShader(fragment);
    }

    ~ShaderSpirv() {
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
