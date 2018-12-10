# PVDiffusionToConductivityTensor

A ParaView plugin to compute conductivity tensors according to the volume constraint method [1]. The plugin is described in the paper **A flexible open-source pipeline for simulating transcranial electric stimulation** by *Benjamin Kalloch, Pierre-Louis Bazin,  Arno Villringer, Bernhard Sehm, and Mario Hlawitschka*.

## How to compile
1. Set up the ParaView environment as described in [2] 
2. Create a folder "build"
3. In the "build" folder execute `cmake ..`
4. Run `make`
5. In ParaView, load the compiled plugin dynamic library (the so-file) under `Tools -> Manage Plugins -> Load New`

###### Useful links
[1] The volume constraint method: https://doi.org/10.1016/j.neuroimage.2005.10.014

[2] Set up the ParaView environment: https://openfoam.org/download/source/third-party-software/