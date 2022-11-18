# Build instructions for the [eX3](https://www.ex3.simula.no/) cluster


In the directory [mfem-build-scripts/ex3](mfem-build-scripts/ex3) you will find
scripts for building mfem-4.5 for the `defq`, `dgx2q`, `mi100q` and `mi210q`.
However, builds for all these partitions are already installed in
`/global/D1/homes/iverh/packages/<partition_name>/mfem<-dbg>-4.5`:

```sh
iverh@srl-login1:~$ tree -L 2 /global/D1/homes/iverh/packages
/global/D1/homes/iverh/packages
├── defq
│   ├── mfem-4.5
│   └── mfem-dbg-4.5
├── dgx2q
│   ├── mfem-4.5
│   └── mfem-dbg-4.5
├── mi100q
│   ├── mfem-4.5
│   └── mfem-dbg-4.5
└── mi210q
    ├── mfem-4.5
    └── mfem-dbg-4.5

12 directories, 0 files

```

To build `cardiac-fibers` for the wanted partition, it should be as simple as running 

```sh
source envsetup-ex3.sh <partition_name>
```

Note that the HIP builds (`mi100q` and `mi210q`) currently have to be built on
one of the nodes that has access to the AMD GPU. For example, to build for
`mi210q`, you can simply do the following:

```sh
salloc -p mi210q --ntasks=1 --gres=rsmi:1
cd <the-cardiac-fibers-directory>
source envsetup-ex3 mi210q
make clean cardiac-fibers
```

### Tested versions

The following setups have been tested on the eX3 cluster:
- A HIP build running on an AMD MI100 GPU on the `mi100q`partition, built with
    * ROCm 5.1.3
    * mfem-4.5 built with HIP support
    * hypre-2.24.0 built with HIP support
    * metis-5.1.0
    * openmpi-4.1.4

- A HIP build running on an AMD MI210 GPU on the `mi210q` partition, built with
    * ROCm 5.1.3
    * mfem-4.5 built with HIP support
    * hypre-2.25.0 built with HIP support
    * metis-5.1.0
    * openmpi-4.1.4

- A CUDA build running on an NVIDIA V100 GPU on the `dgx2q` partition, built with 
    * Cuda Toolkit 10.1.243
    * mfem-4.5 built with CUDA support
    * hypre-2.24.0 built with CUDA support
    * metis-5.1.0
    * openmpi-4.1.4

- A normal CPU build running on the `defq` partition, built with
    * OpenMP support (`CARDIAC_FIBERS_HAS_OPENMP=YES`)
    * mfem-4.5 built with OpenMP support
    * hypre-2.25.0
    * metis-5.1.0
    * openmpi-4.1.4

The builds of mfem used in the tested versions were built using the Slurm
scripts in [mfem-build-scripts/ex3](mfem-build-scripts/ex3). The environment
setup for the tested versions can be found in the script
[envsetup-ex3.sh](envsetup-ex3.sh).
