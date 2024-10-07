
# Direct Volume Rendering with Empty Space Skipping using Partitioned Distance Maps

This is an open-source implementation of [Accelerating Transfer Function Update for Distance Map based Volume Rendering](https://arxiv.org/abs/2407.21552).

Please consider citing this work if you use it:

```
@misc{rauter2024acceleratingtransferfunctionupdate,
      title={Accelerating Transfer Function Update for Distance Map based Volume Rendering}, 
      author={Michael Rauter and Lukas Zimmermann and Markus Zeilinger},
      year={2024},
      eprint={2407.21552},
      archivePrefix={arXiv},
      primaryClass={cs.GR},
      url={https://arxiv.org/abs/2407.21552}, 
}
```

## Features

- GLSL shader for direct volume rendering via raycasting using PDM based ESS
- GLSL compute shaders for computing partitioned occupancy and distance maps
- support for loading volume datasets (SimpleITK supported formats)
- config.ini file in OML format for specifying many different aspects of the demo, e.g.,
    - speficy volume dataset and intensity transfer function (2nd TF can be specified and can be used to blend between both)
    - window size/rendering resolution
    - set windowed or fullscreen mode
    - blocksize parameter
    - number of intensity partitions (= number of PDMs), and special partition for rho_min
    - first person camera mode (WASD+mouse) or demo mode (rotating volume, fixed camera)
    - and more...
- transfer function specification is done in separate .itf file with simple format
- toggle visualization of volume bounding cube, ESS acceleration structure, gradient information

## Limitations
- requires an NVIDIA GPU (tested with RTX 4090, RTX 3060 and GTX 1080 TI), does not work with AMD GPUs (shaders seem to compile fine, but rendered image is black - unfortunately, I do not have an AMD GPU available for investigation)
- transfer function is currently only defined with max 256 control points, TF domain range is stretched to volume intensity range

## Dependencies

- SimpleITK
- Vulkan SDK (used for spirv compilation of GLSL shaders)
- dependencies inside the repository:
    - glad
    - glfw
    - glm
    - toml


## Building


1. build SimpleITK source (https://github.com/SimpleITK/SimpleITK) using SimpleITK's SuperBuild (see instructions here: https://simpleitk.readthedocs.io/en/master/building.html)
1. Vulkan is needed for spirv compilation, get it from here: https://vulkan.lunarg.com/
1. `git clone <this repository>`
1. use CMake or Visual Studio's integrated CMake build tools to build DVR_ESS_PDM

executable will be built into directory `bin/DVR_ESS_PDM/Release` respectively `bin/DVR_ESS_PDM/Debug`
    
## Running

make sure to copy the config file you want to use from the configs folder into the directory of the executable, also make sure that you set config parameters `filename_input_volume` (filepath to volumetric image) and `filetype_input_volume` (SimpleITK input format identifier string).

run using the default `config.ini` file (don't forget to copy it to executable directory):
```bash
DVR_ESS_PDM.exe
```
or
```bash
DVR_ESS_PDM.exe <config_ini_file>
```


### Keyboard Shortcuts
- keys W,A,S,D [only in first person camera mode] ... move camera forward, backward, left, right

- key 1 ... cycle visualization type:
  - volume rendering
  - volume bounding box front face
  - volume bb back face
  - ESS acceleration structure (depending on selected ESS type)
  - voxel gradient

- key 2 ... cycle dithering mode (none, procedurally generated, noise texture [default])

- key 3 ... toggle early ray termination (off, on [default])

- key 4 ... cycle ESS type (none, block ESS, distance ESS [default, PDM method, unless it was deactivated in config-ini file (via `disable_partitioned_distancemap_optimization` option) - in this case distance map will be computed as in paper by Deakin et. al.])

- key 5 ... incorporate light source into volume rendering - how strong should lightsource color be incorporated (light source is not incorporated [default] to max integration in 10 steps)

- key 6 ... cycle gradient type (on the fly computation, precomputed) - initial value set via config-ini file with `precompute_gradient` option, if not set on the fly computation will be default

- key 7 ... enable tricubic texture filtering (off [default], on)

- key 8 ... toggle reusing Normals http://www.marcusbannerman.co.uk/articles/VolumeRendering.html (off [default], on) - has only an effect when light source is incorporated into rendering

- key 9 ... configure diffuse light factor (how much diffuse light is used)  (1.0 [default] in steps of 1.0 up to 10)

- key 0 ... configure specular light factor (how much specular light is used)  (1.0 [default] in steps of 1.0 up to 10)

- key F1 ... if no TF interpolation is just happening and TF1 is not already the active TF, start interpolation to TF1 

- key F2 ... if no TF interpolation is just happening and TF2 is not already the active TF, start interpolation to TF2

### Datasets

volume datasets can be obtained from:

- https://klacansky.com/open-scivis-datasets/
- http://digimorph.org/

### transfer function file format

text encoded file, formatting:
```
applyGammaCorrection:0
gamma:1.21

x_0=0,a_0,r_0,g_0,b_0
x_1,a_1,r_1,g_1,b_1
.
.
.
x_n-1,a_n-1,r_n-1,g_n-1,b_n-1
x_n=255,a_n,r_n,g_n,b_n
```
if applyGammaCorrection is set to 1, gamma correction is applied in shader to transparency contribution of a TF-mapped voxel transparency using gamma value gamma, otherwise gamma is ignored
list of x_n,a_n,r_n,g_n,b_n the control points of the itf and must be unsigned char values, x_0 must be 0, x_n must be 255, any other x_n needs to be in (0,255) and in order and no duplicates!
note: did not have the time to extend this to value ranges > [0,255]
itf is scaled to intensity range of volume

## License
See [LICENSE](LICENSE).
for licenses of third-party dependencies used inside this repository, see files in [lib_extern](https://github.com/CalamityMichL/DVR_ESS_PDM/tree/main/libs_extern) subdirectory licensing.
