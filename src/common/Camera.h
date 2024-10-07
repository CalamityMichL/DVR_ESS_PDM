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
 * code derived from: https://github.com/JoeyDeVries/LearnOpenGL/blob/master/includes/learnopengl/camera.h
 *
 */

#pragma once

#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Camera {
public:
    enum Camera_Movement {
        FORWARD,
        BACKWARD,
        LEFT,
        RIGHT,
        UP,
        DOWN
    };
private:
    glm::vec3 camera_pos;
    glm::vec3 camera_front;
    glm::vec3 camera_up;
    glm::vec3 camera_right;
    glm::vec3 world_up;
    double yaw;
    double pitch;
public:
    Camera() {
        ResetCameraPos();        
    };
    ~Camera() {};

    // returns the view matrix calculated using Euler Angles and the LookAt Matrix
    glm::mat4 GetViewMatrix()
    {
        return glm::lookAt(camera_pos, camera_pos + camera_front, camera_up);
    }

    glm::vec3 GetCameraPos() {
        return camera_pos;
    }

    void SetCameraPos(glm::vec3 pos) {
        camera_pos = pos;
        UpdateCameraVectors();
    }

    void SetCameraPitch(float pitch) {
        this->pitch = pitch;
        UpdateCameraVectors();
    }

    void ResetCameraPos() {
        this->camera_pos = glm::vec3(0.0, 0.0, 2.0);
        this->camera_front = glm::vec3(0.0, 0.0, -1.0);
        this->camera_up = glm::vec3(0.0, 1.0, 0.0);
        this->camera_right = glm::vec3(1.0, 0.0, 0.0);
        this->world_up = glm::vec3(0.0, 1.0, 0.0);
        this->yaw = -90.0f;
        this->pitch = 0.0f;
    }

    void ProcessKeyboard(Camera_Movement direction, float velocity)
    {
        if (direction == FORWARD)
            this->camera_pos += this->camera_front * velocity;
        if (direction == BACKWARD)
            this->camera_pos -= this->camera_front * velocity;
        if (direction == LEFT)
            this->camera_pos -= this->camera_right * velocity;
        if (direction == RIGHT)
            this->camera_pos += this->camera_right * velocity;
        if (direction == UP)
            this->camera_pos += this->camera_up * velocity;
        if (direction == DOWN)
            this->camera_pos -= this->camera_up * velocity;
    }

    void ProcessMouseMovement(double xoffset, double yoffset, bool constrain_pitch = true, float mouse_sensitivity = 0.125f) {
        xoffset *= mouse_sensitivity;
        yoffset *= mouse_sensitivity;

        this->yaw += xoffset;
        this->pitch += yoffset;

        if (constrain_pitch)
        {
            if (this->pitch > 85)
                this->pitch = 85;
            if (this->pitch < -85)
                this->pitch = -85;
        }
        
        UpdateCameraVectors();
    }

private:
    void UpdateCameraVectors()
    {        
        glm::vec3 front;
        front.x = (float)cos(glm::radians(yaw)) * (float)cos(glm::radians(pitch));
        front.y = (float)sin(glm::radians(pitch));
        front.z = (float)sin(glm::radians(yaw)) * (float)cos(glm::radians(pitch));
        camera_front = glm::normalize(front);        
        camera_right = glm::normalize(glm::cross(camera_front, world_up));
        camera_up = glm::normalize(glm::cross(camera_right, camera_front));
    }
};
