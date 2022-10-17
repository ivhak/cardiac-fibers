// Copyright (C) 2022 Iver Håkonsen
//
// ldrb-gpu is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ldrb-gpu is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ldrb-gpu.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include <fstream>
#include <limits>

#include "mfem.hpp"
#include "calculus.hpp"
#include "util.hpp"
#include "ldrb-gpu.hpp"

using namespace mfem;


int main(int argc, char *argv[])
{

    Options opts;

    // Set program defaults
    opts.mesh_file = NULL;
    opts.output_dir = "./out";
    opts.device_config = "hip";
    opts.verbose = 0;
    opts.prescribed_apex = Vector(3);
    opts.paraview = false;

    opts.alpha_endo =  40.0;
    opts.alpha_epi  = -50.0;
    opts.beta_endo  = -65.0;
    opts.beta_epi   =  25.0;

    opts.base_attr = BASE;
    opts.epi_attr  = EPI;
    opts.lv_attr   = LV_ENDO;
    opts.rv_attr   = RV_ENDO;

    // Parse command-line options
    OptionsParser args(argc, argv);
    args.AddOption(&opts.verbose,
            "-v", "--verbose",
            "Be verbose");
    args.AddOption(&opts.mesh_file,
            "-m", "--mesh",
            "Mesh file to use", true);
    args.AddOption(&opts.prescribed_apex,
            "-a", "--apex",
            "Coordinate of apex, space separated list: 'x y z'.", true);
    args.AddOption(&opts.output_dir,
            "-o", "--out",
            "Directory for output files.");
    args.AddOption(&opts.device_config,
            "-d", "--device",
            "Device configuration string, see Device::Configure().");
    args.AddOption(&opts.paraview,
            "-p",  "--paraview",
            "-np", "--no-paraview",
            "Save data files for ParaView (paraview.org) visualization.");
    args.AddOption(&opts.alpha_endo,
            "-ao", "--alpha-endo",
            "Alpha angle in endocardium.");
    args.AddOption(&opts.alpha_epi,
            "-ai", "--alpha-epi",
            "Alpha angle in epicardium.");
    args.AddOption(&opts.beta_endo,
            "-bo", "--beta-endo",
            "Beta angle in endocardium.");
    args.AddOption(&opts.beta_epi,
            "-bi", "--beta-epi",
            "Beta angle in epicardium.");
    args.AddOption(&opts.base_attr,
            "-base", "--base-attribute",
            "Base attribute");
    args.AddOption(&opts.epi_attr,
            "-epi", "--epi-attribute",
            "Epicaridum attribute");
    args.AddOption(&opts.lv_attr,
            "-lv", "--lv-attribute",
            "Left ventricle endocardium attribute");
    args.AddOption(&opts.rv_attr,
            "-rv", "--rv-attribute",
            "Right ventricle endocardium attribute");
    args.Parse();

    if (!args.Good()) {
        args.PrintUsage(std::cout);
        exit(1);
    }

    if (opts.verbose > 1)
        args.PrintOptions(std::cout);

    // Set the basename of the mesh
    opts.mesh_basename = remove_extension(basename(std::string(opts.mesh_file)));

    // Make sure the output direcory exists
    mksubdir(opts.output_dir);

    struct timespec t0, t1;

    // Enable hardware devices such as GPUs, and programming models such as
    // HIP, CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(opts.device_config);
    if (opts.verbose > 1)
        device.Print();

    // Load the mesh
    clock_gettime(CLOCK_MONOTONIC, &t0);
    Mesh mesh(opts.mesh_file, 1, 1);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    if (opts.verbose > 1) {
        std::cout << "Loaded meshfile '" << opts.mesh_file << "' "
             << "consisting of " << mesh.GetNV() << " vertices "
             << "and " << mesh.GetNE() << " elements" << std::endl;
    }

    if (opts.verbose)
        log_timing(std::cout, "Mesh load", timespec_duration(t0, t1));

    // Set the apex node based on the prescribed apex coordinate.
    int apex = 0;
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        apex = find_apex_vertex(mesh, opts.prescribed_apex);
        clock_gettime(CLOCK_MONOTONIC, &t1);

        if (opts.verbose > 1) {
            double *closest_vertex = mesh.GetVertex(apex);
            std::cout << std::setprecision(2)
                 << "Found closest vertex to prescribed apex "
                 << "("  << opts.prescribed_apex[0]
                 << ", " << opts.prescribed_apex[1]
                 << ", " << opts.prescribed_apex[2]
                 << ") at "
                 << "("  << closest_vertex[0]
                 << ", " << closest_vertex[1]
                 << ", " << closest_vertex[2]
                 << ")." << std::endl;

        }
        if (opts.verbose)
            log_timing(std::cout, "Find apex", timespec_duration(t0, t1));
    }

    // Set up the vertex to element table
    Table *vertex_to_element_table = mesh.GetVertexToElementTable();

    int nattr = mesh.bdr_attributes.Max();
    Array<int> ess_bdr(nattr);          // Essential boundaries
    Array<int> nonzero_ess_bdr(nattr);  // Essential boundaries with nonzero value
    Array<int> zero_ess_bdr(nattr);     // Essential boudnaries with zero value

    // Solve the Laplace equation from EPI (1.0) to (opts.rv_attr union opts.rv_attr) (0.0)
    // and calculate the gradients
    GridFunction x_phi_epi;
    double *grad_phi_epi = (double *)malloc(3*mesh.GetNV()*sizeof(double));
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_attr    -1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.epi_attr-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.rv_attr-1] = 1;
        zero_ess_bdr[opts.rv_attr-1] = 1;

        laplace(&x_phi_epi, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, -1, opts.verbose);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_epi", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        calculate_gradients(grad_phi_epi, x_phi_epi, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_epi", timespec_duration(t0, t1));

    }

    // Solve the Laplace equation from opts.rv_attr (1.0) to (opts.rv_attr union EPI) (0.0)
    // and calculate the gradients
    GridFunction x_phi_lv;
    double *grad_phi_lv = (double *)malloc(3*mesh.GetNV()*sizeof(double));
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_attr    -1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.rv_attr-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_attr    -1] = 1;
        zero_ess_bdr[opts.rv_attr-1] = 1;

        laplace(&x_phi_lv, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, -1, opts.verbose);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_lv", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        calculate_gradients(grad_phi_lv, x_phi_lv, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_lv", timespec_duration(t0, t1));
    }


    // Solve the Laplace equation from opts.rv_attr (1.0) to (opts.rv_attr union EPI) (0.0)
    // and calculate the gradients
    GridFunction x_phi_rv;
    double *grad_phi_rv = (double *)malloc(3*mesh.GetNV()*sizeof(double));
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_attr    -1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.rv_attr-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_attr    -1] = 1;
        zero_ess_bdr[opts.rv_attr-1] = 1;

        laplace(&x_phi_rv, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, -1, opts.verbose);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_rv", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        calculate_gradients(grad_phi_rv, x_phi_rv, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_rv", timespec_duration(t0, t1));
    }

    // Solve the Laplace equation from opts.base_attr (1.0) to APEX (0.0)
    // and calculate the gradients
    GridFunction x_psi_ab;
    double *grad_psi_ab = (double *)malloc(3*mesh.GetNV()*sizeof(double));
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        ess_bdr = 0;
        ess_bdr[opts.base_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.base_attr-1] = 1;

        zero_ess_bdr = 0;

        laplace(&x_psi_ab, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, apex, opts.verbose);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "psi_ab", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        calculate_gradients(grad_psi_ab, x_psi_ab, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_psi_ab", timespec_duration(t0, t1));
    }

    delete vertex_to_element_table;


    // Calculate the fiber orientation

    const double *phi_epi = x_phi_epi.Read();
    const double *phi_lv  = x_phi_lv.Read();
    const double *phi_rv  = x_phi_rv.Read();
    const double *psi_ab  = x_psi_ab.Read();

    double *F = (double *)malloc(3 * mesh.GetNV() * sizeof(double));
    double *S = (double *)malloc(3 * mesh.GetNV() * sizeof(double));
    double *T = (double *)malloc(3 * mesh.GetNV() * sizeof(double));

    clock_gettime(CLOCK_MONOTONIC, &t0);
    define_fibers(
            mesh,
            phi_epi,      phi_lv,      phi_rv,      psi_ab,
            grad_phi_epi, grad_phi_lv, grad_phi_rv, grad_psi_ab,
            opts.alpha_endo, opts.alpha_epi, opts.beta_endo, opts.beta_epi,
            F, S, T
    );
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(std::cout, "define_fibers", timespec_duration(t0, t1));


#ifdef DEBUG
    GridFunction grad_phi_epi_gf, grad_phi_lv_gf, grad_phi_rv_gf, grad_psi_ab_gf;
    vertex_vector_to_grid_function(mesh, grad_phi_epi, &grad_phi_epi_gf);
    vertex_vector_to_grid_function(mesh, grad_phi_lv,  &grad_phi_lv_gf);
    vertex_vector_to_grid_function(mesh, grad_phi_rv,  &grad_phi_rv_gf);
    vertex_vector_to_grid_function(mesh, grad_psi_ab,  &grad_psi_ab_gf);
#endif

    GridFunction F_gf, S_gf, T_gf;
    vertex_vector_to_grid_function(mesh, F, &F_gf);
    vertex_vector_to_grid_function(mesh, S, &S_gf);
    vertex_vector_to_grid_function(mesh, T, &T_gf);

    // Save the mesh and solutions
    {

        // Output the normal solutions in the mfem subdirectory
        std::string mfem_output_dir(opts.output_dir);
        mfem_output_dir += "/mfem";
        mksubdir(mfem_output_dir);

#ifdef DEBUG
        std::string debug_dir = mfem_output_dir + "/debug";
        mksubdir(debug_dir);

        debug_print_to_file(grad_phi_lv, 3*mesh.GetNV(), debug_dir, "/grad_phi_lv.txt");
        debug_print_to_file(grad_phi_rv, 3*mesh.GetNV(), debug_dir, "/grad_phi_rv.txt");
        debug_print_to_file(grad_psi_ab, 3*mesh.GetNV(), debug_dir, "/grad_psi_ab.txt");
        debug_print_to_file(F,           3*mesh.GetNV(), debug_dir, "/F.txt");
        debug_print_to_file(S,           3*mesh.GetNV(), debug_dir, "/S.txt");
        debug_print_to_file(T,           3*mesh.GetNV(), debug_dir, "/T.txt");

#endif

        // Save the MFEM mesh
        std::string mesh_out(mfem_output_dir);
        mesh_out += "/" + opts.mesh_basename + ".mesh";
        std::ofstream mesh_ofs(mesh_out.c_str());
        mesh_ofs.precision(8);
        mesh.Print(mesh_ofs);

#ifdef DEBUG
        // Save the solutions
        save_solution(&x_phi_epi, debug_dir, opts.mesh_basename, "_phi_epi.gf");
        save_solution(&x_phi_lv,  debug_dir, opts.mesh_basename, "_phi_lv.gf");
        save_solution(&x_phi_rv,  debug_dir, opts.mesh_basename, "_phi_rv.gf");
        save_solution(&x_psi_ab,  debug_dir, opts.mesh_basename, "_psi_ab.gf");

        save_solution(&grad_phi_epi_gf, debug_dir, opts.mesh_basename, "_grad_phi_epi.gf");
        save_solution(&grad_phi_lv_gf,  debug_dir, opts.mesh_basename, "_grad_phi_lv.gf");
        save_solution(&grad_phi_rv_gf,  debug_dir, opts.mesh_basename, "_grad_phi_rv.gf");
        save_solution(&grad_psi_ab_gf,  debug_dir, opts.mesh_basename, "_grad_psi_ab.gf");
#endif

        save_solution(&F_gf,  mfem_output_dir, opts.mesh_basename, "_F.gf");
        save_solution(&S_gf,  mfem_output_dir, opts.mesh_basename, "_S.gf");
        save_solution(&T_gf,  mfem_output_dir, opts.mesh_basename, "_T.gf");

        // Save in paraview as well
        ParaViewDataCollection *pd = NULL;
        if (opts.paraview) {
            std::string paraview_path(opts.output_dir);
            paraview_path += "/paraview";
            mksubdir(paraview_path);

            pd = new ParaViewDataCollection(opts.mesh_basename, &mesh);
            pd->SetPrefixPath(paraview_path);
#ifdef DEBUG
            pd->RegisterField("grad phi epi", &grad_phi_epi_gf);
            pd->RegisterField("grad phi lv",  &grad_phi_lv_gf);
            pd->RegisterField("grad phi rv",  &grad_phi_rv_gf);
            pd->RegisterField("grad psi ab",  &grad_psi_ab_gf);
            pd->RegisterField("phi epi", &x_phi_epi);
            pd->RegisterField("phi lv",  &x_phi_lv);
            pd->RegisterField("phi rv",  &x_phi_rv);
            pd->RegisterField("psi ab",  &x_psi_ab);
#endif
            pd->RegisterField("F", &F_gf);
            pd->RegisterField("S", &S_gf);
            pd->RegisterField("T", &T_gf);
            pd->SetLevelsOfDetail(1);
            pd->SetDataFormat(VTKFormat::BINARY);
            pd->SetHighOrderOutput(false);
            pd->Save();
        }

    }
}

