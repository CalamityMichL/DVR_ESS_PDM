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

#include <glad/gl.h>

class RenderGeometry {
public:
	GLuint VAO_cube = 0;
	GLuint VBO_cube = 0;
	GLuint EBO_cube = 0;

private:
	float cube_vertices[24] = {
		-0.5f, -0.5f,  0.5f,
		 0.5f, -0.5f,  0.5f,
		 0.5f,  0.5f,  0.5f,
		-0.5f,  0.5f,  0.5f,
		-0.5f, -0.5f, -0.5f,
		 0.5f, -0.5f, -0.5f,
		 0.5f,  0.5f, -0.5f,
		-0.5f,  0.5f, -0.5f
	};

	unsigned int cube_indices[36] = {
		0, 1, 2, 2, 3, 0,
		5, 4, 7, 7, 6, 5,
		4, 5, 1, 1, 0, 4,
		6, 7, 3, 3, 2, 6,
		5, 6, 2, 2, 1, 5,
		7, 4, 0, 0, 3, 7
	};

public:
	RenderGeometry()
	{
		glGenVertexArrays(1, &VAO_cube);
		glBindVertexArray(VAO_cube);

		glGenBuffers(1, &VBO_cube);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_cube);
		glBufferData(GL_ARRAY_BUFFER, 96, cube_vertices, GL_STATIC_DRAW);

		glGenBuffers(1, &EBO_cube);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_cube);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, 144, cube_indices, GL_STATIC_DRAW);
	}

	~RenderGeometry()
	{
		glDeleteVertexArrays(1, &VAO_cube);
		glDeleteBuffers(1, &VBO_cube);
		glDeleteBuffers(1, &EBO_cube);
	}
};
