// Copyright (C) 2022 Iver Håkonsen
//
// cardiac-fibers is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License along with
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

int main(int argc, char **argv)
{
    const char *out_prefix;
    const char *mesh_prefix;
    const char *fiber_prefix;
    const char *trans_fiber_prefix;
    int np = 1;

    OptionsParser args(argc, argv);
    args.AddOption(&out_prefix,         "-o",  "--out",               "Output file", true);
    args.AddOption(&mesh_prefix,        "-m",  "--mesh",              "Input mesh", true);
    args.AddOption(&fiber_prefix,       "-f",  "--fibers",            "Input fiber directions", true);
    args.AddOption(&trans_fiber_prefix, "-t",  "--transverse-fibers", "Input transverse fiber directions", true);
    args.AddOption(&np,                 "-np", "--num-processes",     "Number of processes");
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(std::cout);
        exit(1);
    }

    // Read a parallel Mesh and GridFunctions for fibers and transverse fibers
    // as it is done in ReadParallelMeshAndGridFunction in GLVis
    Array<Mesh *>         mesh_array(np);
    Array<GridFunction *> fibers_array(np);
    Array<GridFunction *> trans_fibers_array(np);

    mesh_array         = NULL;
    fibers_array       = NULL;
    trans_fibers_array = NULL;

    std::ostringstream msg;

    std::cout << "Converting from MFEM to CARP" << std::endl;
    struct timespec start, end, t0, t1;
    util::timing::tick(&start);

    util::timing::tick(&t0);
    if (np > 1) {
        msg << "Reading meshes '"
            << mesh_prefix
            << "{"   << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0
            << "..." << std::setfill('0') <<  std::setw(PAD_DIGITS) << np-1;
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

        mesh_array[p] = new Mesh(mesh_file, 1, 0, 1);
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
    msg.str("");

    util::timing::tick(&t0);
    if (np > 1) {
        msg << "Reading transverse fiber files '"
            << trans_fiber_prefix
            << "{"   << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0
            << "..." << std::setfill('0') <<  std::setw(PAD_DIGITS) << np-1
            << "}'";
    } else {
        msg << "Reading transverse fiber file '"
            << trans_fiber_prefix << "." << std::setfill('0') <<  std::setw(PAD_DIGITS) << 0;
    }

    for (int p = 0; p < np; p++) {
        std::ostringstream trans_fiber_fname;
        trans_fiber_fname << trans_fiber_prefix << '.' << std::setfill('0') << std::setw(PAD_DIGITS) << p;
        ifgzstream trans_fiber_file(trans_fiber_fname.str().c_str());
        if (!trans_fiber_file) {
            std::cerr << "Could not open transverse fibers file "
                << trans_fiber_fname.str() << '!' << std::endl;
            exit(1);
        }
        trans_fibers_array[p] = new GridFunction(mesh_array[p], trans_fiber_file);
    }
    util::timing::tick(&t1);
    util::logging::timestamp(std::cout, msg.str(),
                             util::timing::duration(t0,t1), 1, '+', LEFT_COL_WIDTH);
    msg.str("");

    Mesh *mesh = new Mesh(mesh_array, np);
    GridFunction *fibers = new GridFunction(mesh, fibers_array, np);
    GridFunction *trans_fibers = new GridFunction(mesh, trans_fibers_array, np);

    // Write out the vertices to <out_prefix>.pts
    util::timing::tick(&t0);
    std::ofstream pts_out;
    std::ostringstream pts_out_fname;
    pts_out_fname << out_prefix << ".pts";
    msg << "Writing points to '" << pts_out_fname.str().c_str();
    pts_out.open(pts_out_fname.str().c_str());
    if (pts_out.fail()) {
        std::cerr << "Could not open output file '" << pts_out_fname.str() << "'" << std::endl;
        exit(1);
    }
    Vector coords;
    mesh->GetVertices(coords);
    const int nv = mesh->GetNV();
    const double *vertices = coords.Read();
    pts_out << nv << std::endl;
    for (int i = 0; i < nv; i++) {
        pts_out << std::setprecision(15)
                << 1000*vertices[0*nv + i] << " "
                << 1000*vertices[1*nv + i] << " "
                << 1000*vertices[2*nv + i] << std::endl;
    }
    pts_out.close();
    util::timing::tick(&t1);
    util::logging::timestamp(std::cout, msg.str(),
                            util::timing::duration(t0,t1), 1, '+', LEFT_COL_WIDTH);
    msg.str("");

    // Write out the elements to <out_prefix>.elems
    util::timing::tick(&t0);
    std::ofstream elems_out;
    std::ostringstream elems_out_fname;
    elems_out_fname << out_prefix << ".elem";
    msg << "Writing elements to '" << elems_out_fname.str().c_str();
    elems_out.open(elems_out_fname.str().c_str());
    if (pts_out.fail()) {
        std::cerr << "Could not open output file '" << elems_out_fname.str().c_str() << "'" << std::endl;
        exit(1);
    }
    const int ne = mesh->GetNE();
    elems_out << ne << std::endl;
    for (int i = 0; i < ne; i++) {
        Array<int> vertices;
        mesh->GetElementVertices(i, vertices);
        elems_out << "Tt "
                  << vertices[3] << " "
                  << vertices[0] << " "
                  << vertices[1] << " "
                  << vertices[2] << " 1" << std::endl;
    }
    elems_out.close();
    util::timing::tick(&t1);
    util::logging::timestamp(std::cout, msg.str(),
                             util::timing::duration(t0,t1), 1, '+', LEFT_COL_WIDTH);
    msg.str("");

    // Write out the elemnts to <out_prefix>.elems
    util::timing::tick(&t0);
    std::ofstream fibers_out;
    std::ostringstream fibers_out_fname;
    fibers_out_fname << out_prefix << ".lon";
    msg << "Writing fibers to '" << fibers_out_fname.str().c_str();
    fibers_out.open(fibers_out_fname.str().c_str());
    if (fibers_out.fail()) {
        std::cerr << "Could not open output file '" << fibers_out_fname.str().c_str() << "'" << std::endl;
        exit(1);
    }
    fibers_out << 2 << std::endl;
    const double *fiber_vals = fibers->Read();
    const double *trans_fiber_vals = trans_fibers->Read();
    for (int i = 0; i < ne; i++) {
        fibers_out << std::setprecision(12)
                   << fiber_vals[3*i+0] << " "
                   << fiber_vals[3*i+1] << " "
                   << fiber_vals[3*i+2] << " "
                   << trans_fiber_vals[3*i+0] << " "
                   << trans_fiber_vals[3*i+1] << " "
                   << trans_fiber_vals[3*i+2] << std::endl;
    }
    fibers_out.close();
    util::timing::tick(&t1);
    util::logging::timestamp(std::cout, msg.str(),
                             util::timing::duration(t0,t1), 1, '+', LEFT_COL_WIDTH);
    msg.str("");

    util::timing::tick(&end);
    util::logging::timestamp(std::cout, "Mesh convert done",
                             util::timing::duration(start,end), 1, '=', LEFT_COL_WIDTH);
}
