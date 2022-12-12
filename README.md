# cardiac-fibers

`cardiac-fibers` is an implementation of the Laplace-Dirichlet Rule-Based
algorithm (LDRB) for assigning myocardial fiber orientations, based on the MFEM
finite-element library.

> Bayer, J.D., Blake, R.C., Plank, G. and Trayanova, N.A., 2012.
> A novel rule-based algorithm for assigning myocardial fiber orientation
> to computational heart models. Annals of biomedical engineering, 40(10),
> pp.2243-2254.(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518842/)

Implementations of the LDRB algorithm already exist, see for example
* [ldrb](https://github.com/finsberg/ldrb): An implementation written in
  python, using the [FEniCS](https://fenicsproject.org/) PDE library.
* [life<sup>x</sup>](https://lifex.gitlab.io/): An implementation written in
  C++, using the [deal.II](https://www.dealii.org/) finite element library.
* [cardioid/fibers](https://github.com/LLNL/cardioid/tree/master/fibers): An
  implementation written in C++, which also uses the MFEM library.

The goal of `cardiac-fibers` is first and foremost to implement the LDRB
algorithm with GPU execution in mind, by leveraging the GPU-support
introduced in MFEM version 4.0. 

## Usage

```
Usage: cardiac-fibers [options] ...
Options:
   -h, --help
        Print this help message and exit.
   -v <int>, --verbose <int>, current value: 0
        Set verbosity level:
            1: Print timestamps for each major computation step.
            2: Print more granular timestamps.
            3: Print additional info.
            4: Print BoomerAmg and solver output.
        If the -t (--time-to-file) flag is passed, the output from level 1 and 2
        will be output to <OUT>/time.txt, where <OUT> is set with the -o (--out) flag.
   -m <string>, --mesh <string>, (required)
        Mesh file to use. See https://mfem.org/mesh-formats/ for a list of suppored formats.
   -pm, --par-mesh, -npm, --no-par-mesh, current option: --no-par-mesh
        Mesh is parallel, i.e, one file for each rank.
        Each rank will read the mesh <MESH_FILE>.<rank>, e.g., mesh.000000.
   -a '<double>...', --apex '<double>...', (required)
        Coordinate of apex, space separated list: 'x y z'.
   -o <string>, --out <string>, current value: .
        Directory for output files.
   -t, --time-to-file, -nt, --no-time-to-file, current option: --no-time-to-file
        Output time log to <OUT>/time.txt rather than stdout, where <OUT> is set with the -o (--out) flag.
   -d <string>, --device <string>, current value: cpu
        Device configuration string, see Device::Configure().
   -p, --save-paraview, -np, --no-save-paraview, current option: --save-paraview
        Save data files for ParaView (paraview.org) visualization.
   -s, --save-mfem, -ns, --no-save-mfem, current option: --no-save-mfem
        Save data files in the native MFEM format.
   -sl, --save-laplacians, -nsl, --no-save-laplacians, current option: --save-laplacians
        Save the Laplacians phi_epi, phi_lv and phi_rv
   -sg, --save-gradients, -nsg, --no-save-gradients, current option: --no-save-gradients
        Save the gradients of phi_epi, phi_lv, phi_rv and psi_ab
   -ao <double>, --alpha-endo <double>, current value: 60
        Alpha angle in endocardium.
   -ai <double>, --alpha-epi <double>, current value: -60
        Alpha angle in epicardium.
   -bo <double>, --beta-endo <double>, current value: -60
        Beta angle in endocardium.
   -bi <double>, --beta-epi <double>, current value: 60
        Beta angle in epicardium.
   -base <int>, --base-id <int>, current value: 1
        Id of the base surface.
   -epi <int>, --epi-id <int>, current value: 2
        Id of the epicardium surface.
   -lv <int>, --lv-id <int>, current value: 3
        Id of the left ventricle endocardium surface.
   -rv <int>, --rv-id <int>, current value: 4
        Id of the right ventricle endocardium surface.
        Set to -1 if there is no right ventricle in the geometry, e.g. for a single ventricle geometry.
   -u <int>, --uniform-refinement <int>, current value: 0
        Perform n levels of uniform refinement on the mesh.
   -fpe, --fibers-per-element, -fpv, --fibers-per-vertex, current option: --fibers-per-element
        Calculate fibers in the L2 space (one fiber per element) or H1 (one fiber per vertex).
   -gamg, --gpu-tuned-amg, -ngamg, --no-gpu-tuned-amg, current option: --no-gpu-tuned-amg
        Tune the BoomerAmg preconditioner for (hopefully) better GPU performance.
```

### Example 1: Patient specific bi-ventricle geometry

This example uses the mesh [heart02](mesh/gmsh/heart02.msh), which has been
generated from the data published alongside the following study:

> Martinez-Navarro, Hector, et al. "High arrhythmic risk in antero-septal acute
> myocardial ischemia is explained by increased transmural reentry occurrence."
> Scientific reports 9.1 (2019): 1-12.(https://www.nature.com/articles/s41598-019-53221-2)

The data is published at <https://ora.ox.ac.uk/objects/uuid:951b086c-c4ba-41ef-b967-c2106d87ee06>.

To calculate the fiber orientations we need the id numbers for the base
surface, the epicardium surface, and the left and right ventricle endocardium
surfaces. These can be found in the `PhysicalNames` section of the mesh:

```sh
$ head -n 11 mesh/gmsh/heart02.msh
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
5
2 1 "base"
2 2 "epicardium"
2 3 "left ventricle endocardium"
2 4 "right ventricle endocardium"
3 1 "myocardium"
$EndPhysicalNames
```

Here, the id numbers happen to be the same as the default ones. The second
piece of information we need is where the apex of the epicardium lies. This is
the vertex that is furthest from the base plane. For this mesh, we happen to
know that the apex lies roughly at `[346.35 1233.74 169.79]`.

With the id numbers of the surfaces set (here with the defaults) and the
location of the apex, we can calculate the fiber orientations. We pass the
`--save-paraview` flag to save the solution in a format that can be opened in
the visualization program [ParaView](https://www.paraview.org/). We also pass
the `--verbose` flag with a value of 3 to print some useful information along
the way. It is possible to generate the fibers either on a per element basis
(in $L^2$), or on a per vertex basis (in $H^1$). Here we choose to do it per vertex.


```sh
$ ./cardiac-fibers \
        --mesh mesh/gmsh/heart02.msh \
        --apex '346.35 1233.74 169.79' \
        --fibers-per-vertex \
        --save-paraview \
        --verbose 3
```

We can visualize the generated fibers in ParaView :
1. In ParaView, open the generated file `./paraview/heart02/heart02.pvd`.
2. Click on `'Apply'` to visualize the mesh.
3. With the mesh applied, click on the `'Glyph'` button. 
   1. In the menu that pops up, select `'Orientation Array' > 'F'` and `'Scale Array' > 'No scale array'`.
   2. To increase the number of fibers being rendered, change the number under
      `'Maximum Number of Sample Points'`.
4. Click on  `'Apply'` to visualize the fibers.
5. Select `'phi epi'` on the color menu to color the fibers according their wall depth.

The result should hopefully be something close to this:

![Rendering of the fiber directions on `mesh/gmsh/heart02.msh`](docs/figures/heart02.png)

### Example 2: Idealized left ventricle

It is also possible to generate fibers for a single ventricle. This example
uses the mesh [lv_ellipsoid.msh](mesh/gmsh/lv_ellipsoid.msh), which was
generated using the
[cardiac_geometries](https://github.com/ComputationalPhysiology/cardiac_geometries)
tool.

The `PhysicalNames` section of mesh file shows the id numbers we need:
```sh
$ head -n 15 mesh/gmsh/lv_ellipsoid.msh
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
8
0 1 "ENDOPT"
0 2 "EPIPT"
1 3 "EPIRING"
1 4 "ENDORING"
2 5 "BASE"
2 6 "ENDO"
2 7 "EPI"
3 8 "MYOCARDIUM"
$EndPhysicalNames
$Nodes
```

Inspection of the mesh in [gmsh](https://gmsh.info/) also shows that the apex
of the epicardium is somewhere in the along the $x$ axis, so we set the
prescribed apex coordiantes to `[-100 0 0]`. Note that since we are working
with a single ventricle geometry, we explicitly mark the right ventricle as non
existent by passing the `--rv-id -1` flag.

The mesh itself is quite coarse, so we can refine it using the
`--uniform-refinement` flag. Here we choose to do one level of uniform
refinement, which refines the original mesh from consisting of 770 vertices
2818 elements, to 4909 vertices and 22544 elements.

This time we choose to generate the fibers on a per-element basis rather than
per vertex. To do this we can pass the `--fibers-per-element` (or `-fpe` for
short) flag, to generate fibers in the `L2` space rather than in `H_1`. But,
this is the default behavior, so it can be omitted.

```sh
$ ./cardiac-fibers \
        --mesh mesh/gmsh/lv_ellipsoid.msh \
        --out out/lv_ellipsoid \
        --apex '-100 0 0' \
        --lv-id 6 \
        --base-id 5 \
        --epi-id 7 \
        --rv-id -1 \
        --uniform-refinement 1 \
        --save-paraview \
        --verbose 3

```


![Rendering of the fiber directions on `mesh/gmsh/lv_ellipsoid.msh`](docs/figures/lv_ellipsoid.png)


## Building

* [Build instructions for MacOS](docs/build-darwin.md)

* [Build instructions for the eX3 cluster](docs/build-ex3.md)

`cardiac-fibers` relies minimally on the following libraries:
- [mfem](https://mfem.org), a finite element discretization library.
- [hypre](https://github.com/hypre-space/hypre), a library of high-performance preconditioners and solvers.
- METIS, a family of multilevel partitioning algorithms
- MPI, a message-passing library


The supplied Makefile expects the following variables to be set:

| Variable          | Description                               | Required |
| ----------------- | ------------------------------------------|----------|
| `MFEM_ROOT`       | Root of the built mfem package            | Yes      |
| `MFEM_DBG_ROOT`   | Root of a debug build of the mfem package | No       |
| `HYPRE_INCDIR`    | Location of hypre headers                 | Yes      |
| `HYPRE_LIBDIR`    | Location of hypre library                 | Yes      |
| `METIS_INCDIR`    | Location of METIS headers                 | Yes      |
| `METIS_LIBDIR`    | Location of METIS library                 | Yes      |

A debug build of `cardiac-fibers`, which can be compiled by running `make DEBUG=YES
cardiac-fibers`, requires that the `MFEM_DBG_ROOT` variable is set to the location of
a debug build of mfem. If `MFEM_DBG_ROOT` is not set, `make` will default to
using `MFEM_ROOT`.

See the environment setup script for the [eX3](https://www.ex3.simula.no/)
cluster in [envsetup-ex3.sh](envsetup-ex3.sh) for a working example.

By default (for non-GPU build), the compiler  is set to `mpic++`. This can be
overwritten with the environment variable `MPI_CXX`. Note that `MPI_CXX` has to
point to an actual MPI compiler wrapper, even when compiling for CUDA and HIP,
as it is used to get the compiler flags and link flags for MPI.

### GPU builds

To build with GPU support, `cardiac-fibers` additionally needs a CUDA compiler
(`nvcc`) to run on NVIDIA GPUs, or a HIP compiler (`hipcc`) to run on AMD
GPUs. Note that the GPU builds require that mfem and hypre are also built with
GPU support. See https://github.com/mfem/mfem/blob/master/INSTALL and
https://hypre.readthedocs.io/en/latest/ch-misc.html#building-the-library.
To compile a CUDA-build of `cardiac-fibers`, set the enviroment variable
`CARDIAC_FIBERS_HAS_CUDA=YES`. Similarly, set `CARDIAC_FIBERS_HAS_HIP=YES` for a HIP-build.

By default, the compiler is set to `nvcc` when `CARDIAC_FIBERS_HAS_CUDA=YES`, and
`hipcc` when `CARDIAC_FIBERS_HAS_HIP=YES`. These can be overwritten with the `CUDA_CXX`
and `HIP_CXX` environment variables, respectively.

When building for AMD GPUs it is also possible enable tracing markers in the
code, by running `make HIP_TRACE=YES cardiac-fibers`. This requires the `roctracer` library.

Device execution will be enabled by default when compiling with either CUDA or HIP support.

### OpenMP

If mfem is built with OpenMP (`MFEM_USE_OPENMP=YES`), you can set
`CARDIAC_FIBERS_USE_OPENMP=YES` to build `cardiac-fibers` with OpenMP support as well. To run
with OpenMP, simply set the flag `--device omp`.


## License
`cardiac-fibers` is free software, licensed under the GNU GPL version 3 or (at your
option) any later version. See the file COPYING for copying conditions.
