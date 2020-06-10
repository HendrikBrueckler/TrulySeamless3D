# TrulySeamless3D – Truly Seamless Volumetric Maps

`TrulySeamless3D` is an implementation of the parametrization sanitization procedure described in Sec. 6 of [The 3D Motorcycle Complex for Structured Volume Decomposition \[Brückler et al. 2021\]](https://arxiv.org/abs/2112.05793) (accepted to Eurographics 2022) distributed under GPLv3. It is a 3D generalization of the 2D method from [Exact Constraint Satisfaction for Truly Seamless Parametrization \[Mandad and Campen 2019\]](http://graphics.cs.uos.de/papers/TrulySeamless_EG2019_Mandad.pdf).

If you make use of `TrulySeamless3D` in your scientific work, please cite our paper:

    @article{DBLP:journals/corr/abs-2112-05793,
        author     = {Hendrik Br{\"{u}}ckler and
                     Ojaswi Gupta and
                     Manish Mandad and
                     Marcel Campen},
        title      = {The 3D Motorcycle Complex for Structured Volume Decomposition},
        journal    = {CoRR},
        volume     = {abs/2112.05793},
        year       = {2021},
        url        = {https://arxiv.org/abs/2112.05793},
        eprinttype = {arXiv},
        eprint     = {2112.05793},
    }

***


## What is TrulySeamless3D?

TrulySeamless3D takes as input a tetrahedral mesh with an *almost* seamless parametrization (seamless up to small numerical inaccuracies, as commonly arising in the generation process) and makes it *truly* seamless. In this way, it enables the use of a wide range of almost-seamless parametrizations for applications that depend on the exact satisfaction of seamlessness constraints.


***

### Dependencies
- Eigen (Not included. Must be installed on your system.)
- OpenVolumeMesh (Included as submodule)
- libHexEx (Included as submodule)

### Building
In root directory

    mkdir build
    cd build
    cmake [-DTRULYSEAMLESS_BUILD_CLI=Off] ..
    make

### Usage
A command-line application is included that reads a tetrahedral mesh including an *almost*-seamless parametrization from a file in .hexex-format, as used and documented in [libHexEx](https://www.graphics.rwth-aachen.de/software/libHexEx/).
It outputs a file in the same format, containing a (numerically sanitized) *truly* seamless version of the input.

After building, the CLI app can be found in ```build/Build/bin``` .
Use as follows:

    ./TrulySeamless3D /path/to/input/almost-seamless-map.hexex [--preserve-cut-graph]

```--preserve-cut-graph``` will avoid the cut graph getting altered during sanitization, i.e. non-identity transitions will not be moved onto other facets. Preserving the cut graph may slightly speedup sanitization if the input cut graph is already simple, but may (strongly) deteriorate performance if it is far from simple (e.g. cutting the mesh into many charts, rather than just one). If this flag is not used, internally a relatively simple cut graph will be computed, and transitions be moved onto it.
The output will be saved to ```/path/to/input/almost-seamless-map.hexex-seamless.hexex```.

### API
Call the constructor ```TrulySeamless3D(tetmesh, parameters)```, providing the tetrahedral mesh (of OpenVolumeMesh type) and the (u,v,w) parameters per tet corner (as a map with key tuples (tetID, vertexID)). Then call sanitize(), and finally retrieve the sanitized parametrization map via ```parameter(CellHandle ch, VertexHandle vh)```.
