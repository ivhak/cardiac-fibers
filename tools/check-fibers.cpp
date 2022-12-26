// Copyright (C) 2022 Iver Håkonsen
//
// cardiac-fibers is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include <iostream>
#include <iomanip>

#include "mfem.hpp"
#include "../util.hpp"

using namespace mfem;

#define PAD_DIGITS 6
#define LEFT_COL_WIDTH 72

void read_mesh_and_fiber(const char *mesh_basename, Array<Mesh*> &mesh_array, Array<GridFunction*> &fibers_array, int np, const char *dir)
{

    std::string mesh_prefix   = std::string(dir) + "/mfem/" + mesh_basename + ".mesh";
    std::string fiber_prefix  = std::string(dir) + "/mfem/" + mesh_basename + "_F.gf";

    std::ostringstream msg;

    struct timespec t0, t1;
    util::timing::tick(&t0);
    if (np > 1) {
        msg << "Reading meshes '"
            << mesh_prefix
            << "{"   << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0
            << "..." << std::setfill('0') <<  std::setw(PAD_DIGITS) << np-1 << "}";
    } else {
        msg << "Reading mesh '"
            << mesh_prefix << "." << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0;
    }
    for (int p = 0; p < np; p++) {
        std::ostringstream mesh_fname;
        mesh_fname << mesh_prefix << '.' << std::setfill('0') << std::setw(PAD_DIGITS) << p;
          named_ifgzstream mesh_file(mesh_fname.str().c_str());
        if (!mesh_file) {
            std::cerr << "Could not open mesh file: " << mesh_fname.str() << '!' << std::endl;
            exit(1);
        }

        mesh_array[p] = new Mesh(mesh_file, 1, 1, 1);
    }

    util::timing::tick(&t1);
    util::logging::timestamp(std::cout, msg.str(),
                             util::timing::duration(t0,t1), 1, '+', LEFT_COL_WIDTH);
    msg.str("");

    util::timing::tick(&t0);
    if (np > 1) {
        msg << "Reading fiber files '"
            << fiber_prefix
            << "{"   << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0
            << "..." << std::setfill('0') <<  std::setw(PAD_DIGITS) << np-1
            << "}'";
    } else {
        msg << "Reading fiber file '"
            << fiber_prefix << "."
            << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0;
    }

    for (int p = 0; p < np; p++) {
        std::ostringstream fiber_fname;
        fiber_fname << fiber_prefix << '.' << std::setfill('0') << std::setw(PAD_DIGITS) << p;
        ifgzstream fiber_file(fiber_fname.str().c_str());
        if (!fiber_file) {
            std::cerr << "Could not fibers solution file "
                << fiber_fname.str() << '!' << std::endl;
            exit(1);
        }
        fibers_array[p] = new GridFunction(mesh_array[p], fiber_file);
    }

    util::timing::tick(&t1);
    util::logging::timestamp(std::cout, msg.str(),
                             util::timing::duration(t0,t1), 1, '+', LEFT_COL_WIDTH);
}

void save_fiber(
    mfem::GridFunction *x,
    std::string const& dir,
    std::string const& base_name,
    std::string const& suffix)
{
    std::string filename(dir);
    filename += "/";
    filename += base_name;
    filename += suffix;
    std::ofstream x_ofs(filename.c_str());
    x_ofs.precision(12);
    x->Save(x_ofs);
}

int main(int argc, char **argv)
{

    const char *mesh_basename;
    const char *serial_dir;
    const char *parallel_dir;
    const char *out;
    int np = 1;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_basename, "-m",  "--mesh-basename", "Basename of the mesh", true);
    args.AddOption(&serial_dir,    "-s",  "--serial",        "Directory hold the serial output", true);
    args.AddOption(&parallel_dir,  "-p",  "--parallel",      "Directory holding the parallel output", true);
    args.AddOption(&out,           "-o",  "--out",           "Output directory", true);
    args.AddOption(&np,            "-np", "--num-processes", "Number of processes in the parallel output", true);
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(std::cout);
        exit(1);
    }

    Mesh *smesh = NULL;
    GridFunction *sfiber = NULL;
    Array<Mesh *> smesh_array(1);
    Array<GridFunction *> sfiber_array(1);

    Array<Mesh *> pmesh_array(np);
    Array<GridFunction *> pfiber_array(np);

    std::cout << "Loading meshes and fibers" << std::endl;
    read_mesh_and_fiber(mesh_basename, smesh_array, sfiber_array, 1,  serial_dir);
    read_mesh_and_fiber(mesh_basename, pmesh_array, pfiber_array, np, parallel_dir);

    smesh  = new Mesh(smesh_array, 1);
    sfiber = new GridFunction(smesh, sfiber_array, 1);

    const int num_elements = smesh->GetNE();

    FiniteElementSpace *s_fes = sfiber->FESpace();

    std::cout << "Reading partitioning" << std::endl;
    int *partitioning = (int *)malloc(num_elements * sizeof(int));
    std::string partitioning_file(parallel_dir);
    partitioning_file += "/mfem/partitioning.txt";
    std::fstream partitioning_in(partitioning_file);
    for (int i = 0; i < num_elements; i++) {
        partitioning_in >> partitioning[i];
    }

    const double **local_fibers = (const double **)malloc(np * sizeof(double *));
    for (int i = 0; i < np; i++) {
        local_fibers[i] = pfiber_array[i]->Read();
    }

    std::cout << "Converting parallel fibers to serial" << std::endl;
    GridFunction pfiber_serial(s_fes);
    int *rank_idx = (int *)calloc(np, sizeof(int));
    double *pfiber_serial_vals = pfiber_serial.Write();
    for (int i = 0; i < num_elements; i++) {
        int element_in_rank = partitioning[i];
        int local_idx = rank_idx[element_in_rank]++;
        pfiber_serial_vals[3*i+0] = local_fibers[element_in_rank][3*local_idx+0];
        pfiber_serial_vals[3*i+1] = local_fibers[element_in_rank][3*local_idx+1];
        pfiber_serial_vals[3*i+2] = local_fibers[element_in_rank][3*local_idx+2];
    }

    // Calculate the "diff" between the serial and parallel solution bu
    // subtracting the parallel from the serial. If correct, the resulting
    // vector should have a magnitude of 0.
    std::cout << "Calculating diff against serial" << std::endl;
    GridFunction diff(s_fes);
    diff = *sfiber;
    diff -= pfiber_serial;

    std::cout << "" << std::endl;
    ParaViewDataCollection *pd = NULL;
    std::string paraview_path(out);
    paraview_path += "/paraview";
    util::fs::mksubdir(paraview_path);

    {
        std::string msg = "Saving ParaView files to '" + paraview_path + "'";
        std::cout << msg << std::endl;
    }
    pd = new ParaViewDataCollection(mesh_basename, smesh);
    pd->SetPrefixPath(paraview_path);
    pd->RegisterField("F_serial",   sfiber);
    pd->RegisterField("F_parallel", &pfiber_serial);
    pd->RegisterField("F_diff",     &diff);
    pd->SetLevelsOfDetail(1);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(false);
    pd->Save();

    std::string mfem_dir(out);
    mfem_dir += "/mfem";
    util::fs::mksubdir(mfem_dir);
    {
        std::string msg = "Saving MFEM files to '" + mfem_dir + "'";
        std::cout << msg << std::endl;
    }
    save_fiber(sfiber,         mfem_dir, mesh_basename, "_F_serial.gf");
    save_fiber(&pfiber_serial, mfem_dir, mesh_basename, "_F_parallel.gf");
    save_fiber(&diff,          mfem_dir, mesh_basename, "_F_diff.gf");

    return 0;
}
