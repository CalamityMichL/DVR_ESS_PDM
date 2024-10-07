# Copyright Michael Rauter
#
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 the "License";
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

set(CMAKE_PROJECT_NAME "DVR_ESS_PDM")

if (CMAKE_BUILD_TYPE)
	set(OUTPUTDIR ${CMAKE_SOURCE_DIR}/bin/${CMAKE_PROJECT_NAME}/${CMAKE_BUILD_TYPE})
	set(EXECUTABLE_OUTPUT_PATH ${OUTPUTDIR})
else()
	if (CMAKE_GENERATOR MATCHES "Visual Studio")
		set(OUTPUTDIR ${CMAKE_SOURCE_DIR}/bin/${CMAKE_PROJECT_NAME}/$(Configuration))
		set(EXECUTABLE_OUTPUT_PATH ${OUTPUTDIR}/../)
	endif()
endif()

set(CMAKE_CURRENT_BINARY_DIR ${OUTPUTDIR})

find_package(SimpleITK REQUIRED)
find_package(Vulkan COMPONENTS glslc)
find_program(glslc_executable NAMES glslc HINTS Vulkan::glslc)

function(compile_shader target)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "ENV;INPUTPATH;INPUTFORMAT;OUTPUTPATH;OUTPUTFORMAT;OUTPUTFILEPREFIX" "SOURCES")
    foreach(source ${arg_SOURCES})
        add_custom_command(            
            OUTPUT ${arg_OUTPUTPATH}/${arg_OUTPUTFILEPREFIX}${source}.${arg_OUTPUTFORMAT}
            DEPENDS ${arg_INPUTPATH}${source}.${arg_INPUTFORMAT}
            DEPFILE ${arg_INPUTPATH}${source}.${arg_INPUTFORMAT}.d
            COMMAND
                ${glslc_executable}
                $<$<BOOL:${arg_ENV}>:--target-env=${arg_ENV}>
                #$<$<BOOL:${arg_OUTPUTFORMAT}>:-mfmt=${arg_OUTPUTFORMAT}>
                #-MD -MF ${arg_INPUTPATH}${source}.${arg_INPUTFORMAT}.d                
                -o ${arg_OUTPUTPATH}/${arg_OUTPUTFILEPREFIX}${source}.${arg_OUTPUTFORMAT}
                ${CMAKE_CURRENT_SOURCE_DIR}/${arg_INPUTPATH}${source}.${arg_INPUTFORMAT}
        )
        target_sources(${target} PRIVATE ${arg_OUTPUTPATH}/${arg_OUTPUTFILEPREFIX}${source}.${arg_OUTPUTFORMAT})
    endforeach()
endfunction()

function(compile_shader_variant target)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "ENV;INPUTPATH;INPUTFORMAT;OUTPUTPATH;OUTPUTFORMAT;OUTPUTFILEPREFIX;OUTPUTFILEPOSTFIX;PREPROCESSOR_DEFINES" "SOURCE_FILE")
    if ("${arg_PREPROCESSOR_DEFINES}" STREQUAL "")
        set(ADDITIONAL_COMMAND_LINE_PARAMS)
    else()
        set(ADDITIONAL_COMMAND_LINE_PARAMS)        
        foreach(elem IN LISTS ${arg_PREPROCESSOR_DEFINES})
            set(ADDITIONAL_COMMAND_LINE_PARAMS ${ADDITIONAL_COMMAND_LINE_PARAMS} -D${elem})
        endforeach()
        #set(ADDITIONAL_COMMAND_LINE_PARAMS "-D\"${arg_PREPROCESSOR_DEFINES}\"")
    endif()
    add_custom_command(
        OUTPUT ${arg_OUTPUTPATH}/${arg_OUTPUTFILEPREFIX}${arg_SOURCE_FILE}${arg_OUTPUTFILEPOSTFIX}.${arg_OUTPUTFORMAT}
        DEPENDS ${arg_INPUTPATH}${arg_SOURCE_FILE}.${arg_INPUTFORMAT}
        DEPFILE ${arg_INPUTPATH}${arg_SOURCE_FILE}.${arg_INPUTFORMAT}.d
        COMMAND
            ${glslc_executable}
            $<$<BOOL:${arg_ENV}>:--target-env=${arg_ENV}>
            #$<$<BOOL:${arg_OUTPUTFORMAT}>:-mfmt=${arg_OUTPUTFORMAT}>
            #-MD -MF ${arg_INPUTPATH}${source}.${arg_INPUTFORMAT}.d                
            -o ${arg_OUTPUTPATH}/${arg_OUTPUTFILEPREFIX}${arg_SOURCE_FILE}${arg_OUTPUTFILEPOSTFIX}.${arg_OUTPUTFORMAT}
            ${ADDITIONAL_COMMAND_LINE_PARAMS}
            ${CMAKE_CURRENT_SOURCE_DIR}/${arg_INPUTPATH}${arg_SOURCE_FILE}.${arg_INPUTFORMAT}
    )
    target_sources(${target} PRIVATE ${arg_OUTPUTPATH}/${arg_OUTPUTFILEPREFIX}${arg_SOURCE_FILE}${arg_OUTPUTFILEPOSTFIX}.${arg_OUTPUTFORMAT})
endfunction()

include_directories("../../libs_extern/glfw/include")
include_directories("../../libs_extern/glad/include")
include_directories("../../libs_extern/glm")

add_executable (DVR_ESS_PDM "DVR_ESS_PDM.cpp" "../../libs_extern/glad/include/KHR/khrplatform.h" "../../libs_extern/glad/include/glad/gl.h" "../../libs_extern/glad/src/gl.c" "../common/Camera.h" "../common/TextureHelper.h" "../LaunchComputeShaders.h" "../LaunchComputeShaders.cpp" "../ITF.h" "../Volume.h" "../common/ShaderSpirv.h" "../common/ComputeShaderSpirv.h")

set(COMPILE_SHADERS TRUE)

if(${COMPILE_SHADERS})
set(preprocessor_defines "VOLUME_DATATYPE_INT16" "VOLUME_DATATYPE_OUPUT_UNORM16")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX 16_from_dtype_int16
    SOURCE_FILE convert_volume_to_unorm
)
set(preprocessor_defines "VOLUME_DATATYPE_UINT16" "VOLUME_DATATYPE_OUPUT_UNORM16")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX 16_from_dtype_uint16
    SOURCE_FILE convert_volume_to_unorm
)
set(preprocessor_defines "VOLUME_DATATYPE_UINT8" "VOLUME_DATATYPE_OUPUT_UNORM16")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX 16_from_dtype_uint8
    SOURCE_FILE convert_volume_to_unorm
)
set(preprocessor_defines "VOLUME_DATATYPE_INT16" "VOLUME_DATATYPE_OUPUT_UNORM8")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX 8_from_dtype_int16
    SOURCE_FILE convert_volume_to_unorm
)
set(preprocessor_defines "VOLUME_DATATYPE_UINT16" "VOLUME_DATATYPE_OUPUT_UNORM8")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX 8_from_dtype_uint16
    SOURCE_FILE convert_volume_to_unorm
)
set(preprocessor_defines "VOLUME_DATATYPE_UINT8" "VOLUME_DATATYPE_OUPUT_UNORM8")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX 8_from_dtype_uint8
    SOURCE_FILE convert_volume_to_unorm
)

set(preprocessor_defines "VOLUME_DATATYPE_UNORM8")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX _dtype_unorm8
    SOURCE_FILE compute_gradient_map
)
set(preprocessor_defines "VOLUME_DATATYPE_UNORM16")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX _dtype_unorm16
    SOURCE_FILE compute_gradient_map
)

set(preprocessor_defines "VOLUME_DATATYPE_UNORM8")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX _dtype_unorm8
    SOURCE_FILE compute_occupancy_map
)
set(preprocessor_defines "VOLUME_DATATYPE_UNORM16")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX _dtype_unorm16
    SOURCE_FILE compute_occupancy_map
)

set(preprocessor_defines "VOLUME_DATATYPE_UNORM8")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX _dtype_unorm8
    SOURCE_FILE compute_POM
)
set(preprocessor_defines "VOLUME_DATATYPE_UNORM16")
compile_shader_variant(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    PREPROCESSOR_DEFINES preprocessor_defines
    OUTPUTFILEPOSTFIX _dtype_unorm16
    SOURCE_FILE compute_POM
)

compile_shader(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT glsl
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX
    SOURCES
        vs_render_volume
        fs_render_volume
)

compile_shader(DVR_ESS_PDM
    ENV opengl
    INPUTPATH ../../shaders/
    INPUTFORMAT comp
    OUTPUTPATH ${OUTPUTDIR}/shaders/
    OUTPUTFORMAT spv
    OUTPUTFILEPREFIX cs_
    SOURCES        
        compute_distance_map
        combine_PDMs
)
endif()

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET DVR_ESS_PDM PROPERTY CXX_STANDARD 20)
endif()

if (CMAKE_BUILD_TYPE)
	set( GLFW_LIB ${CMAKE_SOURCE_DIR}/libs_extern/glfw/lib/${CMAKE_BUILD_TYPE}/glfw3.lib )
else()
	if (CMAKE_GENERATOR MATCHES "Visual Studio")
		set( GLFW_LIB ${CMAKE_SOURCE_DIR}/libs_extern/glfw/lib/$(Configuration)/glfw3.lib )
	endif()	
endif()

target_link_libraries(DVR_ESS_PDM ${GLFW_LIB} ${SimpleITK_LIBRARIES})

add_custom_command(TARGET DVR_ESS_PDM POST_BUILD    
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/itfs/" "${CMAKE_CURRENT_BINARY_DIR}/itfs/"
)

