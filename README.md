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
        Mesh file to use.
        See https://mfem.org/mesh-formats/ for a list of suppored formats.
   -a '<double>...', --apex '<double>...', (required)
        Coordinate of apex, space separated list: 'x y z'.
   -o <string>, --out <string>, current value: ./out
        Directory for output files.
   -t, --time-to-file, -nt, --no-time-to-file, current option: --no-time-to-file
        Output time log to <OUT>/time.txt rather than stdout,
        where <OUT> is set with the -o (--out) flag.
   -d <string>, --device <string>, current value: hip
        Device configuration string, see Device::Configure().
   -p, --save-paraview, -np, --no-save-paraview, current option: --no-save-paraview
        Save data files for ParaView (paraview.org) visualization.
   -s, --save-mfem, -ns, --no-save-mfem, current option: --save-mfem
        Save data files in the native MFEM format.
   -it <double>, --interpolation-tolerance <double>, current value: 0.1
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
        Id of the base surface
   -epi <int>, --epi-id <int>, current value: 2
        Id of the epicardium surface
   -lv <int>, --lv-id <int>, current value: 3
        Id of the left ventricle endocardium surface.
   -rv <int>, --rv-id <int>, current value: 4
        Id of the right ventricle endocardium surface.
        Set to -1 if there is no right ventricle in the geometry,
        e.g. for a single ventricle geometry.
   -s <int>, --solver <int>, current value: 0
        Solver to use. Options are:
            0: HyprePCG (See mfem::HyprePCG)
            1: CGSolver (See mfem::CGSolver)
   -gamg, --gpu-tuned-amg, -ngamg, --no-gpu-tuned-amg, current option: --no-gpu-tuned-amg
        Tune the BoomerAmg preconditioner for (hopefully) better GPU performance.
```

## License
`ldrb-gpu` is free software, licensed under the GNU GPL version 3 or (at your option) any later version. See the file COPYING for copying conditions.
