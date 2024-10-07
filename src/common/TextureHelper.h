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
 */

#pragma once

#include <glad/gl.h>
#include <glm/glm.hpp>

class TextureHelper {
public:
    static GLuint* CreateTexture3D(GLint internalformat, GLsizei width, GLsizei height, GLsizei depth, GLenum format, GLenum type, const void* dataArray, GLint filterMode, GLint wrappingMode, GLuint* textureIDs, int numberOfTexturesToCreate = 1)
    {
        GLfloat borderColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

        glGenTextures(numberOfTexturesToCreate, textureIDs);
        for (int i = 0; i < numberOfTexturesToCreate; i++)
        {
            glBindTexture(GL_TEXTURE_3D, textureIDs[i]);
            glTexImage3D(GL_TEXTURE_3D, 0, internalformat, width, height, depth, 0, format, type, dataArray);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, wrappingMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, wrappingMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, wrappingMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, filterMode);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, filterMode);
            glTexParameterfv(GL_TEXTURE_3D, GL_TEXTURE_BORDER_COLOR, borderColor);
                
        }
        return textureIDs;
    }

    static GLuint* CreateTexture2D(GLint internalformat, GLsizei width, GLsizei height, GLenum format, GLenum type, void* dataArray, GLint filterMode, GLint wrappingMode, GLuint* textureIDs, int numberOfTexturesToCreate = 1)
    {
        GLfloat borderColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

        glGenTextures(numberOfTexturesToCreate, textureIDs);
        for (int i = 0; i < numberOfTexturesToCreate; i++)
        {
            glBindTexture(GL_TEXTURE_2D, textureIDs[i]);
            glTexImage2D(GL_TEXTURE_2D, 0, internalformat, width, height, 0, format, type, dataArray);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrappingMode);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrappingMode);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filterMode);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filterMode);
            glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

        }
        return textureIDs;
    }

    static GLuint* CreateTexture1D(GLint internalformat, GLsizei size, GLenum format, GLenum type, void* dataArray, GLint filterMode, GLint wrappingMode, GLuint* textureIDs, int numberOfTexturesToCreate = 1)
    {
        GLfloat borderColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

        glGenTextures(numberOfTexturesToCreate, textureIDs);
        for (int i = 0; i < numberOfTexturesToCreate; i++)
        {
            glBindTexture(GL_TEXTURE_1D, textureIDs[i]);
            glTexImage1D(GL_TEXTURE_1D, 0, internalformat, size, 0, format, type, dataArray);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, wrappingMode);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, wrappingMode);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, filterMode);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, filterMode);
            glTexParameterfv(GL_TEXTURE_1D, GL_TEXTURE_BORDER_COLOR, borderColor);

        }
        return textureIDs;
    }

    static void DeleteTexture(GLuint textureID)
    {
        if (textureID == NULL)
            return;

        glDeleteTextures(1, &textureID);
        textureID = NULL;
    }
    static void DeleteTextures(GLsizei n, GLuint* textureIDs)
    {
        if (textureIDs == NULL)
            return;

        glDeleteTextures(n, textureIDs);

        for (GLsizei i=0; i < n; i++)
            textureIDs[i] = NULL;
    }

};
