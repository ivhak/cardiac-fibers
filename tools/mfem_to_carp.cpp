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

using namespace mfem;

#define PAD_DIGITS 6

int main(int argc, char **argv)
{
    const char *out_prefix = "";
    const char *mesh_prefix = "";
    const char *fiber_prefix = "";
    const char *trans_fiber_prefix = "";
    int np;

    OptionsParser args(argc, argv);
    args.AddOption(&out_prefix,         "-o",  "--out",               "Output file");
    args.AddOption(&mesh_prefix,        "-m",  "--mesh",              "Input mesh");
    args.AddOption(&fiber_prefix,       "-f",  "--fibers",            "Input fiber directions");
    args.AddOption(&trans_fiber_prefix, "-t",  "--transverse-fibers", "Input transverse fiber directions");
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

    int read_err = 0;
    for (int p = 0; p < np; p++) {
        std::ostringstream mesh_fname;
        mesh_fname << mesh_prefix << '.' << std::setfill('0') << std::setw(PAD_DIGITS) << p;
          named_ifgzstream mesh_file(mesh_fname.str().c_str());
        if (!mesh_file) {
            std::cerr << "Could not open mesh file: " << mesh_fname.str() << '!' << std::endl;
            read_err = 1;
            break;
        }

        mesh_array[p] = new Mesh(mesh_file, 1, 0, 1);


        std::ostringstream fiber_fname;
        fiber_fname << fiber_prefix << '.' << std::setfill('0') << std::setw(PAD_DIGITS) << p;
        ifgzstream fiber_file(fiber_fname.str().c_str());
        if (!fiber_file) {
            std::cerr << "Could not fibers solution file "
                << fiber_fname.str() << '!' << std::endl;
            read_err = 2;
            break;
        }
        fibers_array[p] = new GridFunction(mesh_array[p], fiber_file);

        std::ostringstream trans_fiber_fname;
        trans_fiber_fname << trans_fiber_prefix << '.' << std::setfill('0') << std::setw(PAD_DIGITS) << p;
        ifgzstream trans_fiber_file(trans_fiber_fname.str().c_str());
        if (!trans_fiber_file) {
            std::cerr << "Could not open transverse fibers file "
                << fiber_fname.str() << '!' << std::endl;
            read_err = 2;
            break;
        }
        trans_fibers_array[p] = new GridFunction(mesh_array[p], trans_fiber_file);
    }

    Mesh *mesh = new Mesh(mesh_array, np);
    GridFunction *fibers = new GridFunction(mesh, fibers_array, np);
    GridFunction *trans_fibers = new GridFunction(mesh, trans_fibers_array, np);

    // Write out the vertices to <out_prefix>.pts
    std::ofstream pts_out;
    std::ostringstream pts_out_fname;
    pts_out_fname << out_prefix << ".pts";
    std::cout << "Writing points to '" << pts_out_fname.str().c_str() << "'" << std::endl;
    pts_out.open(pts_out_fname.str().c_str());
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

    // Write out the elements to <out_prefix>.elems
    std::ofstream elems_out;
    std::ostringstream elems_out_fname;
    elems_out_fname << out_prefix << ".elem";
    std::cout << "Writing elements to '" << elems_out_fname.str().c_str() << "'" << std::endl;
    elems_out.open(elems_out_fname.str().c_str());
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

    // Write out the elemnts to <out_prefix>.elems
    std::ofstream fibers_out;
    std::ostringstream fibers_out_fname;
    fibers_out_fname << out_prefix << ".lon";
    std::cout << "Writing fibers to '" << fibers_out_fname.str().c_str() << "'" << std::endl;
    fibers_out.open(fibers_out_fname.str().c_str());
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





}
