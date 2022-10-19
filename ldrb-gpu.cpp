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

void laplace(
    GridFunction *x,
    Mesh& mesh,
    Array<int> &ess_bdr,
    Array<int> &nonzero_ess_bdr,
    Array<int> &zero_ess_bdr,
    int apex,
    int verbose)
{

    // Determine the list of true (i.e. parallel conforming) essential boundary
    // dofs, defined by the boundary attributes marked as essential (Dirichlet)
    // and converting to a list of true dofs..
    Array<int> ess_tdof_list;
    MFEM_ASSERT(mesh.bdr_attributes.Size() != 0, "Boundary size cannot be zero.");

    FiniteElementSpace *fespace = x->FESpace();

    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Set up the linear form b(.) which corresponds to the right-hand
    // side of the FEM linear system, which in this case is (0,phi_i) where
    // phi_i are the basis functions in fespace.
    LinearForm b(fespace);
    ConstantCoefficient zero(0.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));
    b.Assemble();

    // Initialize x with initial guess of zero, which
    // satisfies the boundary conditions.
    *x = 0.0;

    // Project the constant value 1.0 to all the essential boundaries marked as nonzero.
    ConstantCoefficient nonzero_bdr(1.0);
    x->ProjectBdrCoefficient(nonzero_bdr, nonzero_ess_bdr);

    // Project the constant value 0.0 to all the essential boundaries marked as zero.
    ConstantCoefficient zero_bdr(0.0);
    x->ProjectBdrCoefficient(zero_bdr, zero_ess_bdr);

    // For the laplacians involving the apex we need to treat the boundary
    // conditions a little different, as it is to be enforced on a single node
    // rather than a whole boundary surface. To do so, we make sure that the
    // apex is en essential true dof, and then we project the wanted value, in
    // this case 0.0, to only that node.
    if (apex >= 0) {
        // Initialize the internal data needed in the finite element space
        fespace->BuildDofToArrays();

        // Make sure the apex is in the list of essential true Dofs
        ess_tdof_list.Append(apex);

        Array<int> node_disp(1);
        node_disp[0] = apex;
        Vector node_disp_value(1);
        node_disp_value[0] = 0.0;

        VectorConstantCoefficient node_disp_value_coeff(node_disp_value);

        x->ProjectCoefficient(node_disp_value_coeff, node_disp);
    }

    // Set up the parallel bilinear form a(.,.) on the finite element space
    // corresponding to the Laplacian operator -Delta, by adding the
    // Diffusion domain integrator.
    BilinearForm a(fespace);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    // Assemble the parallel bilinear form and the corresponding linear system.
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, *x, b, A, X, B);

    // Solve the linear system A X = B.
    // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, verbose > 1 ? 1 : 0, 1000, 1e-15, 0.0);

    // Recover the grid function corresponding to X.
    a.RecoverFEMSolution(X, b, *x);
}

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

    opts.geom_has_rv = true;

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
            "Right ventricle endocardium attribute. "
            "Set to -1 if there is no right ventricle in the geometry, "
            "e.g. for a single ventricle geometry.");
    args.Parse();

    if (!args.Good()) {
        args.PrintUsage(std::cout);
        exit(1);
    }

    if (opts.verbose > 1)
        args.PrintOptions(std::cout);

    if (opts.rv_attr == -1) {
        opts.geom_has_rv = false;
    }

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

    H1_FECollection fec(1, mesh.Dimension());

    // Set up two finite element spaces: one for scalar values (1D) , and one for 3d vectors
    FiniteElementSpace fespace1d(&mesh, &fec);
    FiniteElementSpace fespace3d(&mesh, &fec, 3, mfem::Ordering::byVDIM);

    // Set up the vertex to element table
    Table *vertex_to_element_table = mesh.GetVertexToElementTable();

    int nattr = mesh.bdr_attributes.Max();
    Array<int> ess_bdr(nattr);          // Essential boundaries
    Array<int> nonzero_ess_bdr(nattr);  // Essential boundaries with nonzero value
    Array<int> zero_ess_bdr(nattr);     // Essential boudnaries with zero value

    // Solve the Laplace equation from EPI (1.0) to (opts.rv_attr union opts.rv_attr) (0.0)
    // and calculate the gradients
    GridFunction x_phi_epi(&fespace1d);
    GridFunction grad_phi_epi(&fespace3d);
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_attr-1] = 1;
        ess_bdr[opts.lv_attr -1] = 1;
        if (opts.geom_has_rv)
            ess_bdr[opts.rv_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.epi_attr-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.lv_attr-1] = 1;
        if (opts.geom_has_rv)
            zero_ess_bdr[opts.rv_attr-1] = 1;

        laplace(&x_phi_epi, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, -1, opts.verbose);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_epi", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        double *grad_phi_epi_vals = grad_phi_epi.Write();
        calculate_gradients(grad_phi_epi_vals, x_phi_epi, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_epi", timespec_duration(t0, t1));

    }

    // Solve the Laplace equation from opts.rv_attr (1.0) to (opts.rv_attr union EPI) (0.0)
    // and calculate the gradients
    GridFunction x_phi_lv(&fespace1d);
    GridFunction grad_phi_lv(&fespace3d);
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_attr-1] = 1;
        ess_bdr[opts.lv_attr -1] = 1;
        if (opts.geom_has_rv)
            ess_bdr[opts.rv_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.lv_attr-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_attr-1] = 1;
        if (opts.geom_has_rv)
            zero_ess_bdr[opts.rv_attr-1] = 1;

        laplace(&x_phi_lv, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, -1, opts.verbose);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_lv", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        double *grad_phi_lv_vals = grad_phi_lv.Write();
        calculate_gradients(grad_phi_lv_vals, x_phi_lv, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_lv", timespec_duration(t0, t1));
    }


    // Solve the Laplace equation from opts.rv_attr (1.0) to (opts.rv_attr union EPI) (0.0)
    // and calculate the gradients
    GridFunction x_phi_rv(&fespace1d);
    GridFunction grad_phi_rv(&fespace3d);
    if (opts.geom_has_rv) {
        clock_gettime(CLOCK_MONOTONIC, &t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_attr    -1] = 1;
        ess_bdr[opts.lv_attr-1] = 1;
        ess_bdr[opts.rv_attr-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.rv_attr-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_attr-1] = 1;
        zero_ess_bdr[opts.lv_attr -1] = 1;

        laplace(&x_phi_rv, mesh, ess_bdr, nonzero_ess_bdr, zero_ess_bdr, -1, opts.verbose);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_rv", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        double *grad_phi_rv_vals = grad_phi_rv.Write();
        calculate_gradients(grad_phi_rv_vals, x_phi_rv, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_rv", timespec_duration(t0, t1));
    } else {
        // If the goemetry does not have a right ventricle, just make a dummy,
        // zero filled GridFunction. The rest of the algorithm _should_ be
        // implementeded in such a way that it still produces the right result.
        x_phi_rv = 0.0;
        grad_phi_rv = 0.0;
    }

    // Solve the Laplace equation from opts.base_attr (1.0) to APEX (0.0)
    // and calculate the gradients
    GridFunction x_psi_ab(&fespace1d);
    GridFunction grad_psi_ab(&fespace3d);
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
        double *grad_psi_ab_vals = grad_psi_ab.Write();
        calculate_gradients(grad_psi_ab_vals, x_psi_ab, mesh, vertex_to_element_table);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_psi_ab", timespec_duration(t0, t1));
    }

    delete vertex_to_element_table;



    // Get read-only pointers to the internal arrays of the laplacian solutions
    const double *phi_epi = x_phi_epi.Read();
    const double *phi_lv  = x_phi_lv.Read();
    const double *phi_rv  = x_phi_rv.Read();
    const double *psi_ab  = x_psi_ab.Read();

    // Get read-only pointers to the internal arrays of the gradients
    const double *grad_phi_epi_vals = grad_phi_epi.Read();
    const double *grad_phi_lv_vals  = grad_phi_lv.Read();
    const double *grad_phi_rv_vals  = grad_phi_rv.Read();
    const double *grad_psi_ab_vals  = grad_psi_ab.Read();

    // Setup GridFunctions to store the fibre directions in
    GridFunction F(&fespace3d);
    GridFunction S(&fespace3d);
    GridFunction T(&fespace3d);

    double *F_vals = F.Write();
    double *S_vals = S.Write();
    double *T_vals = T.Write();


    // Calculate the fiber orientation
    clock_gettime(CLOCK_MONOTONIC, &t0);
    define_fibers(
            mesh,
            phi_epi,      phi_lv,      phi_rv,      psi_ab,
            grad_phi_epi_vals, grad_phi_lv_vals, grad_phi_rv_vals, grad_psi_ab_vals,
            opts.alpha_endo, opts.alpha_epi, opts.beta_endo, opts.beta_epi,
            F_vals, S_vals, T_vals
    );
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(std::cout, "define_fibers", timespec_duration(t0, t1));

    // Save the mesh and solutions
    {

        // Output the normal solutions in the mfem subdirectory
        std::string mfem_output_dir(opts.output_dir);
        mfem_output_dir += "/mfem";
        mksubdir(mfem_output_dir);

#ifdef DEBUG
        std::string debug_dir = mfem_output_dir + "/debug";
        mksubdir(debug_dir);

        debug_print_to_file(grad_phi_epi_vals, 3*mesh.GetNV(), debug_dir, "/grad_phi_epi.txt");
        debug_print_to_file(grad_phi_lv_vals,  3*mesh.GetNV(), debug_dir, "/grad_phi_lv.txt");
        debug_print_to_file(grad_phi_rv_vals,  3*mesh.GetNV(), debug_dir, "/grad_phi_rv.txt");
        debug_print_to_file(grad_psi_ab_vals,  3*mesh.GetNV(), debug_dir, "/grad_psi_ab.txt");
        debug_print_to_file(F,                 3*mesh.GetNV(), debug_dir, "/F.txt");
        debug_print_to_file(S,                 3*mesh.GetNV(), debug_dir, "/S.txt");
        debug_print_to_file(T,                 3*mesh.GetNV(), debug_dir, "/T.txt");

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

        save_solution(&grad_phi_epi, debug_dir, opts.mesh_basename, "_grad_phi_epi.gf");
        save_solution(&grad_phi_lv,  debug_dir, opts.mesh_basename, "_grad_phi_lv.gf");
        save_solution(&grad_phi_rv,  debug_dir, opts.mesh_basename, "_grad_phi_rv.gf");
        save_solution(&grad_psi_ab,  debug_dir, opts.mesh_basename, "_grad_psi_ab.gf");
#endif

        save_solution(&F,  mfem_output_dir, opts.mesh_basename, "_F.gf");
        save_solution(&S,  mfem_output_dir, opts.mesh_basename, "_S.gf");
        save_solution(&T,  mfem_output_dir, opts.mesh_basename, "_T.gf");

        // Save in paraview as well
        ParaViewDataCollection *pd = NULL;
        if (opts.paraview) {
            std::string paraview_path(opts.output_dir);
            paraview_path += "/paraview";
            mksubdir(paraview_path);

            pd = new ParaViewDataCollection(opts.mesh_basename, &mesh);
            pd->SetPrefixPath(paraview_path);
#ifdef DEBUG
            pd->RegisterField("grad phi epi", &grad_phi_epi);
            pd->RegisterField("grad phi lv",  &grad_phi_lv);
            pd->RegisterField("grad phi rv",  &grad_phi_rv);
            pd->RegisterField("grad psi ab",  &grad_psi_ab);
            pd->RegisterField("phi epi", &x_phi_epi);
            pd->RegisterField("phi lv",  &x_phi_lv);
            pd->RegisterField("phi rv",  &x_phi_rv);
            pd->RegisterField("psi ab",  &x_psi_ab);
#endif
            pd->RegisterField("F", &F);
            pd->RegisterField("S", &S);
            pd->RegisterField("T", &T);
            pd->SetLevelsOfDetail(1);
            pd->SetDataFormat(VTKFormat::BINARY);
            pd->SetHighOrderOutput(false);
            pd->Save();
        }

    }
}

