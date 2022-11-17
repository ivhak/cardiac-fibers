# A GPU-based implementation of the Laplace-Dirichlet Rule-Based (LDRB) algorithm for assigning myocardial fiber orientations

A GPU-based implementation of the LDRB algorithm using the MFEM library (https://mfem.org/).

> Bayer, J.D., Blake, R.C., Plank, G. and Trayanova, N.A., 2012.
> A novel rule-based algorithm for assigning myocardial fiber orientation
>to computational heart models. Annals of biomedical engineering, 40(10),
pp.2243-2254.(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518842/)

## Usage

```
Usage: ldrb-gpu [options] ...
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
   -a '<double>...', --apex '<double>...', (required)
        Coordinate of apex, space separated list: 'x y z'.
   -o <string>, --out <string>, current value: ./out
        Directory for output files.
   -t, --time-to-file, -nt, --no-time-to-file, current option: --no-time-to-file
        Output time log to <OUT>/time.txt rather than stdout, where <OUT> is set with the -o (--out) flag.
   -d <string>, --device <string>, current value: cpu
        Device configuration string, see Device::Configure().
   -p, --save-paraview, -np, --no-save-paraview, current option: --no-save-paraview
        Save data files for ParaView (paraview.org) visualization.
   -s, --save-mfem, -ns, --no-save-mfem, current option: --save-mfem
        Save data files in the native MFEM format.
   -it <double>, --interpolation-tolerance <double>, current value: 1e-12
        Tolerance for LDRB interpolations.
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
   -dg, --discontinuous-galerkin, -ndg, --no-discontinuous-galerkin, current option: --no-discontinuous-galerkin
        Calculate fibers in the DG0 space (one fiber per element) rather than H1 (one fiber per vertex).
   -gamg, --gpu-tuned-amg, -ngamg, --no-gpu-tuned-amg, current option: --no-gpu-tuned-amg
        Tune the BoomerAmg preconditioner for (hopefully) better GPU performance.
```

### Example: Idealized left ventricle

This example uses the mesh [lv_ellipsoid.msh](mesh/gmsh/lv_ellipsoid.msh),
which was generated using the
[cardiac_geometries](https://github.com/ComputationalPhysiology/cardiac_geometries)
tool. The `PhysicalNames` section of mesh file shows the id numbers we need:
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
prescribed apex coordiantes to `[-100 0 0 ]`. Note that since we are working
with a single ventricle geometry, we explicitly mark the right ventricle as non
existent by passing the `--rv-id -1` flag. The mesh itself is quite coarse, so
we can refine it using the `--uniform-refinement` flag. Here we choose to do
two levels of uniform refinement, which refines the original mesh from
consisting of 770 vertices 2818 elements, to consisting of 34569 vertices and
180352 elements. With all the information we need, we can finally run the following command

```sh
$ ./ldrb-gpu --mesh mesh/gmsh/lv_ellipsoid.msh \
             --apex '-100 0 0' \
             --lv-id 6 \
             --base-id 5 \
             --epi-id 7 \
             --rv-id -1 \
             --uniform-refinement 2 \
             --out out/lv_ellipsoid \
             --save-paraview 

```

Opening the resulting solution
`./out/lv_ellipsoid/paraview/lv_ellipsoid/lv_ellipsoid.pvd` in paraview gives
the following results:

![Rendering of the fiber directions on `mesh/gmsh/lv_ellipsoid.msh`](docs/figures/lv_ellipsoid.png)


## Building

`ldrb-gpu` relies minimally on the following libraries:
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
| `MPI_INCDIR`      | Location of MPI headers                   | Yes      |
| `MPI_LIBDIR`      | Location of MPI library                   | Yes      |

A debug build of `ldrb-gpu`, which can be compiled by running `make DEBUG=YES
ldrb-gpu`, requires that the `MFEM_DBG_ROOT` variable is set to the location of
a debug build of mfem. If you do not have a debug build of mfem, this can be
circumvented by setting it to the same as `MFEM_ROOT`, for example by running
`MFEM_DBG_ROOT=$MFEM_ROOT make DEBUG=YES ldrb-gpu`.

See the environment setup script for the [eX3](https://www.ex3.simula.no/)
cluster in [envsetup-ex3.sh](envsetup-ex3.sh) for a working example.

By default (for non-GPU build), the compiler  is set to `mpiCC`. This can be
overwritten with the environment variable `MPI_CXX`.

### GPU builds

To build with GPU support, `ldrb-gpu` additionally needs a CUDA compiler
(`nvcc`) to run on NVIDIA GPUs, or a HIP compiler (`hipcc`) to run on AMD
GPUs. Note that the GPU builds require that mfem and hypre are also built with
GPU support. See https://github.com/mfem/mfem/blob/master/INSTALL and
https://hypre.readthedocs.io/en/latest/ch-misc.html#building-the-library.
To compile a CUDA-build of `ldrb-gpu`, set the enviroment variable
`LDRB_HAS_CUDA=YES`. Similarly, set `LDRB_HAS_HIP=YES` for a HIP-build.

By default, the compiler is set to `nvcc` when `LDRB_HAS_CUDA=YES`, and
`hipcc` when `LDRB_HAS_HIP=YES`. These can be overwritten with the `CUDA_CXX`
and `HIP_CXX` environment variables, respectively.

When building for AMD GPUs it is also possible enable tracing markers in the
code, by running `make HIP_TRACE=YES ldrb-gpu`. This requires the `roctracer` library.

Device execution will be enabled by default when compiling with either CUDA or HIP support.

### OpenMP

If mfem is built with OpenMP (`MFEM_USE_OPENMP=YES`), you can set
`LDRB_USE_OPENMP=YES` to build `ldrb-gpu` with OpenMP support as well. To run
with OpenMP, simply set the flag `--device omp`.

### Tested versions

The following setups have been tested on the [eX3](https://www.ex3.simula.no/) cluster:
- A HIP build running on an AMD MI100 GPU on the `mi100q`partition, built with
    * ROCm 5.1.3
    * mfem-4.5 built with HIP support
    * hypre-2.24.0 built with HIP support
    * metis-5.1.0
    * openmpi-4.1.4

- A CUDA build running on an NVIDIA V100 GPU on the `dgx2q` partition, built with 
    * Cuda Toolkit 10.1.243
    * mfem-4.5 built with CUDA support
    * hypre-2.24.0 built with CUDA support
    * metis-5.1.0
    * openmpi-4.1.4

- A normal CPU build running on the `defq` partition, built with
    * OpenMP support (`LDRB_HAS_OPENMP=YES`)
    * mfem-4.5 built with OpenMP support
    * hypre-2.25.0
    * metis-5.1.0
    * openmpi-4.1.4

The builds of mfem used in the tested versions were built using the Slurm
scripts in [mfem-build-scripts/ex3](mfem-build-scripts/ex3). The environment
setup for the tested versions can be found in the script
[envsetup-ex3.sh](envsetup-ex3.sh).

## License
`ldrb-gpu` is free software, licensed under the GNU GPL version 3 or (at your
option) any later version. See the file COPYING for copying conditions.
