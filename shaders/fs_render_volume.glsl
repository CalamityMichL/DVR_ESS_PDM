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
 * code derived from: https://github.com/LDeakin/VkVolume/blob/master/shaders/volume_render.frag
 *
 */

#version 460

#pragma shader_stage(fragment)

precision highp int;
precision highp float;

layout (location = 1) in vec4 position;
layout (location = 0) out vec4 out_color;

layout (binding = 0) uniform sampler3D texVolume;

layout (binding = 1) uniform sampler3D texGradient;
layout (binding = 2) uniform sampler2D texNoise;
layout (binding = 3) uniform mediump sampler2D texITF;
layout (binding = 4) uniform mediump usampler3D texOccupancyMap;
layout (binding = 5) uniform mediump usampler3D texDistanceMap;

layout(std140, binding = 1) uniform UBO {
    mat4 model;                     // aligned offset 0
    mat4 view;                      // aligned offset 64
    mat4 inverseModelView;          // aligned offset 128
    float sampling_factor;          // aligned offset 192
    int renderContentType;          // aligned offset 196
    int ditheringMode;              // aligned offset 200
    bool earlyRayTermination;       // aligned offset 204
    int essType;                    // aligned offset 208
    int gradientType;               // aligned offset 212
    bool itf_gamma_correction;      // aligned offset 216
    float gamma;                    // aligned offset 220
    bool enable_bicubic_filtering;  // aligned offset 224
    bool reuse_normals;             // aligned offset 228
    float lightshading_factor;      // aligned offset 232
    vec3 lightDirection;            // aligned offset 240 // must be memory aligned at 16 * x -> 16*15 == 240
    float light_ambient_factor;     // aligned offset 252
    float light_diffuse_factor;     // aligned offset 256
    float light_specular_factor;    // aligned offset 260
    float light_specular_shininess; // aligned offset 264
} ubo;

#define MAX_SEGMENTS 64

// Ray
struct Ray {
    vec3 origin;
    vec3 direction;
};

// Axis-aligned bounding box
struct AABB {
    vec3 max;
    vec3 min;
};

// Slab method for ray-box intersection
void ray_box_intersection(Ray ray, AABB box, out float t_0, out float t_1)
{
    vec3 direction_inv = 1.0 / ray.direction;
    vec3 t_top = direction_inv * (box.max - ray.origin);
    vec3 t_bottom = direction_inv * (box.min - ray.origin);
    vec3 t_min = min(t_top, t_bottom);
    vec2 t = max(t_min.xx, t_min.yz);
    t_0 = max(0.0, max(t.x, t.y));
    vec3 t_max = max(t_top, t_bottom);
    t = min(t_max.xx, t_max.yz);
    t_1 = min(t.x, t.y);
}

float getVolumeIntensityValue(vec3 samplingPos, vec3 volumeSize, vec3 texelWidth, bool forceNoTriCubic)
{
    if (!ubo.enable_bicubic_filtering || forceNoTriCubic)
    {
        float voxelValue = texture(texVolume, samplingPos.xyz).r;
        return voxelValue;
    }
    else
    {
        vec3 coord_hg = samplingPos * volumeSize - vec3(0.5f, 0.5f, 0.5);
        vec3 a = coord_hg - floor(coord_hg);

        vec3 w_0 = (1.0/6.0)*(a*(a*(-a + 3.0) - 3.0) + 1.0);
        vec3 w_1 = (1.0/6.0)*(a*a*(3.0*a - 6.0) + 4.0);
        vec3 w_2 = (1.0/6.0)*(a*(a*(-3.0*a + 3.0) + 3.0) + 1.0);
        vec3 w_3 = (1.0/6.0)*(a*a*a);

        vec3 g_0 = w_0 + w_1;
        vec3 g_1 = w_2 + w_3;
        vec3 h_0 = 1 - (w_1/(w_0 + w_1)) + a;
        vec3 h_1 = 1 + (w_3/(w_2 + w_3)) - a;

        vec3 texCoord000 = samplingPos.xyz + vec3(-texelWidth.x * h_0.x, -texelWidth.y * h_0.y, -texelWidth.z * h_0.z);
        vec3 texCoord100 = samplingPos.xyz + vec3(+texelWidth.x * h_1.x, -texelWidth.y * h_0.y, -texelWidth.z * h_0.z);
        vec3 texCoord010 = samplingPos.xyz + vec3(-texelWidth.x * h_0.x, +texelWidth.y * h_1.y, -texelWidth.z * h_0.z);
        vec3 texCoord110 = samplingPos.xyz + vec3(+texelWidth.x * h_1.x, +texelWidth.y * h_1.y, -texelWidth.z * h_0.z);
        vec3 texCoord001 = samplingPos.xyz + vec3(-texelWidth.x * h_0.x, -texelWidth.y * h_0.y, +texelWidth.z * h_1.z);
        vec3 texCoord101 = samplingPos.xyz + vec3(+texelWidth.x * h_1.x, -texelWidth.y * h_0.y, +texelWidth.z * h_1.z);
        vec3 texCoord011 = samplingPos.xyz + vec3(-texelWidth.x * h_0.x, +texelWidth.y * h_1.y, +texelWidth.z * h_1.z);
        vec3 texCoord111 = samplingPos.xyz + vec3(+texelWidth.x * h_1.x, +texelWidth.y * h_1.y, +texelWidth.z * h_1.z);


        float voxelValue000 = texture(texVolume, texCoord000).r;
        float voxelValue100 = texture(texVolume, texCoord100).r;
        float voxelValue010 = texture(texVolume, texCoord010).r;
        float voxelValue110 = texture(texVolume, texCoord110).r;
        float voxelValue001 = texture(texVolume, texCoord001).r;
        float voxelValue101 = texture(texVolume, texCoord101).r;
        float voxelValue011 = texture(texVolume, texCoord011).r;
        float voxelValue111 = texture(texVolume, texCoord111).r;

        g_0 = 1.0f - g_0;
        float l1 = mix(voxelValue000, voxelValue100, g_0.x);
        float l2 = mix(voxelValue010, voxelValue110, g_0.x);
        float l3 = mix(voxelValue001, voxelValue101, g_0.x);
        float l4 = mix(voxelValue011, voxelValue111, g_0.x);
        float l5 = mix(l1,l2,g_0.y);
        float l6 = mix(l3,l4,g_0.y);
        float l7 = mix(l5,l6,g_0.z);

        return l7;
    }
}

vec4 getGradient(vec3 samplingPos, vec3 volumeSize, float voxelValue, int gradientType)
{
    vec4 gradient;
    if (gradientType == 0) // on-the-fly approximation of gradient using tetrahedron technique
    {
        vec2 k = ivec2(-1.0f, 1.0f);
        vec3 gradientVec = 0.25*(
        k.xyy * getVolumeIntensityValue(samplingPos + (k.xyy + 0.5)/vec3(volumeSize), volumeSize, 1.0/volumeSize, true) +
        k.yyx * getVolumeIntensityValue(samplingPos + (k.yyx + 0.5)/vec3(volumeSize), volumeSize, 1.0/volumeSize, true) +
        k.yxy * getVolumeIntensityValue(samplingPos + (k.yxy + 0.5)/vec3(volumeSize), volumeSize, 1.0/volumeSize, true) +
        k.xxx * getVolumeIntensityValue(samplingPos + (k.xxx + 0.5)/vec3(volumeSize), volumeSize, 1.0/volumeSize, true));
        gradient = vec4(normalize(gradientVec), length(gradientVec));        
    }    
    else if (gradientType == 1) // use precomputed gradient
    {
        gradient = vec4(normalize(texture(texGradient, samplingPos + 0.5/vec3(volumeSize)).xyz), texture(texGradient, samplingPos + 0.5/vec3(volumeSize)).a);
    }

    return gradient;
}

void main()
{
    mat4 model = ubo.model;
    mat4 view = ubo.view;
    mat4 inverseModelView = ubo.inverseModelView;

    float sampling_factor = ubo.sampling_factor;

    int renderContentType = ubo.renderContentType;
    int ditheringMode = ubo.ditheringMode;
    bool earlyRayTermination = ubo.earlyRayTermination;
    int essType = ubo.essType;
    float lightshading_factor = ubo.lightshading_factor;    
    int gradientType = ubo.gradientType;
    bool reuse_normals = ubo.reuse_normals;
    bool itf_gamma_correction = ubo.itf_gamma_correction;
    float gamma = ubo.gamma;
    vec3 lightDirection = ubo.lightDirection;
    float light_ambient_factor = ubo.light_ambient_factor;
    float light_diffuse_factor = ubo.light_diffuse_factor;
    float light_specular_factor = ubo.light_specular_factor;
    float light_specular_shininess = ubo.light_specular_shininess;

    ivec2 dimITF = textureSize(texITF, 0).xy;

    vec4 camPosInViewCoords = vec4(view[0][3],view[1][3],view[2][3],1.0f);
    vec4 camPosInObjectCoords = inverseModelView * camPosInViewCoords;
    vec3 ray_origin = camPosInObjectCoords.xyz/camPosInObjectCoords.w;

    vec3 bb_bf_positions = position.xyz;
    vec3 ray_origin_to_bf = bb_bf_positions - ray_origin;
    vec3 ray_direction = normalize(ray_origin_to_bf);

    float t_0, t_1;
    vec3 maxBB = vec3(0.5f);
    vec3 minBB = vec3(-0.5f);
    Ray casting_ray = Ray(ray_origin, ray_direction);
    AABB bounding_box = AABB(maxBB, minBB);
    ray_box_intersection(casting_ray, bounding_box, t_0, t_1);

    vec3 ray_start = (ray_origin + ray_direction * t_0 - minBB) / (maxBB - minBB);
    vec3 ray_stop = (ray_origin + ray_direction * t_1 - minBB) / (maxBB - minBB);

    vec3 ray = ray_stop - ray_start;
    float ray_length = length(ray);
    vec3 ray_dir = normalize(ray);

    if (renderContentType == 1)
    {
        out_color = vec4(ray_start, 1.0f); return;
    }
    else if (renderContentType == 2)
    {
        out_color = vec4(ray_stop, 1.0f); return;
    }

    // termination of rays (pixels) where the ray is degenerated (e.g. start and endpoint coinciding, ...) - this will prevent performance regressions
    if (isinf(ray_length) || isnan(ray_length) || ray_length > 1.732051f) // sqrt(3)
        return;

    ivec3 volumeSize = textureSize(texVolume, 0).xyz;
    int dim_max = max(max(volumeSize.x, volumeSize.y), volumeSize.z);

    int n_steps = int(ceil(float(dim_max) * ray_length * sampling_factor));
    vec3 step_vector = normalize(ray_dir) * ray_length / (float(n_steps) - 1.0f);
    float sampling_factor_inv = 1.0f / sampling_factor;

    // ray "dithering" mode which adds small random substep (either derived from frag coords or from noise texture)
    float random = 0.0f;
    if (ditheringMode == 1)
        random = fract(sin(gl_FragCoord.x * 12.9898 + gl_FragCoord.y * 78.233) * 43758.5453);
    else if (ditheringMode == 2)
        random = texture(texNoise, gl_FragCoord.xy / 256).r;

    vec3 samplingPos;

    // Empty space skipping initialization
    ivec3 dim_distance_map = textureSize(texDistanceMap, 0).xyz;
    ivec3 dim_distance_map_1 = dim_distance_map - 1;
    ivec3 block_size = ivec3(ceil(vec3(volumeSize) / vec3(dim_distance_map)));
    vec3 volume_to_distance_map = vec3(volumeSize) / (vec3(block_size) * vec3(dim_distance_map));
    vec3 step_dist_texel = step_vector * vec3(volumeSize) / vec3(block_size);
    vec3 step_dist_texel_inv = 1.0f / step_dist_texel;
    bool skipping = true;
    int i_min = 0; // current furthest "occupied" step

    out_color = vec4(0.0f);

    vec4 color = vec4(0.0f);
    vec4 prevGrad = vec4(0.0f);

    int distance_map_idx = (ray_dir.z < 0 ? 1 : 0) + (ray_dir.y < 0 ? 2 : 0) + (ray_dir.x < 0 ? 4 : 0);

    vec3 volumeSizeF = vec3(volumeSize);

    for (int i = 0; i < n_steps;)
    {
        vec3 samplingPos = ray_start + float(i + random) * step_vector;        

        if (essType > 0 && skipping && renderContentType != 3 ) {
            vec3 pos_distance_map = (volume_to_distance_map * samplingPos);
            vec3 u = pos_distance_map * vec3(dim_distance_map);
            ivec3 u_i = ivec3(clamp(floor(u), ivec3(0), dim_distance_map_1));
            uint dist;
            if (essType == 1)
                dist = uint(texelFetch(texOccupancyMap, u_i, 0).x);
            else if (essType == 2)
                dist = uint(texelFetch(texDistanceMap, u_i, 0).x);
            vec3 r = clamp(u_i - u, -1.0, 0.0);

            if (dist > 0u) {
                vec3 i_delta_xyz = vec3(0);
                if (essType == 1) // Skip with "block empty space skipping"
                    i_delta_xyz = (step(0.0f, step_dist_texel_inv) + r) * step_dist_texel_inv;
                else if (essType == 2) // Skip with "chebyshev empty space skipping"
                    i_delta_xyz = (step(0.0f, -step_dist_texel_inv) + sign(step_dist_texel_inv) * float(dist) + r) * step_dist_texel_inv;

                int i_delta = max(1, int(ceil(min(min(i_delta_xyz.x, i_delta_xyz.y), i_delta_xyz.z))));
                i += i_delta;
            }
            else
            {
                int i_delta_backwards = -int(ceil(sampling_factor));
                skipping = false;
                i = max(i + i_delta_backwards, i_min);
            }
        }
        else
        {
            float voxelValue = getVolumeIntensityValue(samplingPos.xyz+0.5/volumeSize, volumeSize, 1.0/volumeSize, false);

            vec4 color = texture(texITF, vec2(voxelValue+0.5/dimITF[0], 0.5/dimITF[1])).rgba;

            if (itf_gamma_correction)
                color.a = pow(color.a, gamma);

            if (lightshading_factor > 0.0)
            {
                vec4 gradient = vec4(0.0f);

                gradient = getGradient(samplingPos, volumeSize, voxelValue, gradientType);

                if (reuse_normals)
                {
                    //reusing Normals http://www.marcusbannerman.co.uk/articles/VolumeRendering.html
                    if (gradient.a < prevGrad.a)
                        gradient.rgb = prevGrad.rgb;
                    prevGrad = gradient;
                }

                vec3 normal = - gradient.xyz;

                vec3 light = (inverseModelView * view * vec4(lightDirection,1.0)).xyz;

                float diffuse = (light_diffuse_factor) * max(0.0, min(1.0, dot(normal.rgb, normalize(light))));

                vec3 viewDir    = normalize(ray_origin - samplingPos);
                vec3 halfwayDir = normalize(light + viewDir);
                float spec = light_specular_factor * pow(max(dot(normal.rgb, halfwayDir), 0.0), light_specular_shininess);

                color.rgb = (1.0f - lightshading_factor) * color.rgb + color.rgb * lightshading_factor * (spec + diffuse + light_ambient_factor);

                if (renderContentType == 4)
                {
                    color = vec4(normal, gradient.a);
                }
            }

            if (renderContentType == 3) // render occupancy map, distance map
            {
                float dist = 0.0f;
                        
                vec3 pos_distance_map = (volume_to_distance_map * samplingPos);
                vec3 u = pos_distance_map * vec3(dim_distance_map);
                ivec3 u_i = ivec3(clamp(floor(u), ivec3(0), dim_distance_map_1));                
                if (essType == 1)
                    dist = texelFetch(texOccupancyMap, u_i, 0).x;
                else if (essType == 2)
                    dist = texelFetch(texDistanceMap, u_i, 0).x;
                    
                color = vec4(vec3(dist), dist/dim_max);
            }
            
            color.a = clamp(1.0f * (1.0f - pow(1.0f - color.a, sampling_factor_inv)), 0.0f, 1.0f);

            color.xyz *= color.a;

            out_color = out_color + (1.0f - out_color.a) * color;

            // Early ray termination
            if (out_color.a > 0.98f && earlyRayTermination) {              
              out_color.a = 1.0f;
              break;
            }

            i_min = i;

            // Move the ray forward
            ++i;
        }
    }
}
