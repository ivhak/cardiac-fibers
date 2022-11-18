# Build instructions for MacOS

> Disclaimer: This has only been tested on MacOS Catalina 10.15.7, and it
> assumes that the user is using the [Homebrew](https://brew.sh/) package
> manager.


## Install dependencies

The dependencies can either be built from source, or we can just install them using Homebrew.

```sh
$ brew install open-mpi hypre metis
```

## Build MFEM

First, download mfem from <https://mfem.org/>. This tutorial will assume
version 4.5. Then we can unpack the source, and step into the directory:

```sh
tar -xzvf mfem-4.5.tgz
cd mfem-4.5

```

Then we can build mfem, either by copying and running the  [build
script](../mfem-build-scripts/build_mfem_darwin.sh) which will install mfem
into `${HOME}/.local/packages/mfem-4.5`, or by running the following commands:

```sh

export MFEM_BUILD_DIR=build-darwin                       # For example
export MFEM_INSTALL_DIR=${HOME}/.local/packages/mfem-4.5 # For example

make BUILD_DIR=${MFEM_BUILD_DIR} config \
    MFEM_USE_MPI=YES \
    MPICXX=mpic++ \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L/usr/local/lib -lHYPRE" \
    HYPRE_OPT="-L/usr/local/include" \
    METIS_LIB="-L/usr/local/lib -lmetis" \
    METIS_OPT="-L/usr/local/include"

# Build mfem, using 4 processes to speed things up. Tailor to your system.
make BUILD_DIR=${MFEM_BUILD_DIR} -j 4

# Quick check that everything works. Note that this may fail if your system
# can't handle `mpirun -np 4`. If this command errors for that reason, just skip
# along.
make BUILD_DIR=${MFEM_BUILD_DIR} check

make BUILD_DIR=${MFEM_BUILD_DIR} install PREFIX=${MFEM_INSTALL_DIR}

```

## Build cardiac-fibers

With mfem installed, we can finally build `cardiac-fibers`. To do so, we need
to set up the required environment variables. If you followed this guide, you
can (hopefully) just run the following command from the root of the
`cardiac-fibers` repo:

```sh
source envsetup-darwin.sh
```

If you deviated from the guide at any point, you can look at the script and
change things accordingly. Then, it should be as simple as running

```sh
make cardiac-fibers
```


