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

#include <string>
#include <unistd.h>
#include <sys/stat.h>

#include "util.hpp"

double timespec_duration(
    struct timespec t0,
    struct timespec t1)
{
    return (t1.tv_sec - t0.tv_sec) +
        (t1.tv_nsec - t0.tv_nsec) * 1e-9;
}

void log_timing(
    std::ostream& out,
    const char *log_string,
    double seconds)
{
    out << "[" << std::left << std::setw(15) << log_string << "]: "
        << std::right << std::fixed << std::setw(12)
        << std::setprecision(6)<< seconds << " s" << std::endl;
}

std::string basename(std::string const& path)
{
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string remove_extension(std::string const& filename)
{
    size_t index_of_extension = filename.find_last_of(".");
    return filename.substr(0, index_of_extension);
}

void mksubdir(std::string const& subdir)
{
    const char *cwd = getcwd(NULL, 0);
    std::string path(cwd);
    path += "/";
    path += subdir;
    mkdir(path.c_str(), 0777);
}

// Save a solution (in form of a GridFunction) to a file named "<prefix><suffix>".
void save_solution(
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
    x_ofs.precision(8);
    x->Save(x_ofs);
}

// Save a solution (in form of a GridFunction) to a file named "<prefix><suffix>".
void save_solution(
    mfem::ParGridFunction *x,
    std::string const& dir,
    std::string const& base_name,
    std::string const& suffix,
    int rank)
{
    std::ostringstream filename;
    filename << dir << "/"
             << base_name
             << suffix << "."
             << std::setfill('0') << std::setw(6) << rank;
    std::ofstream x_ofs(filename.str().c_str());
    x_ofs.precision(8);
    x->Save(x_ofs);
}

// Find the vertex closest to the prescribed apex, Euclidean distance.
int find_apex_vertex(mfem::Mesh& mesh, mfem::Vector& apex)
{
    int apex_vertex = 0;
    double distance = std::numeric_limits<double>::max();
    for (int i = 0; i < mesh.GetNV(); i++) {
        double *vertices = mesh.GetVertex(i);
        double this_distance = (vertices[0]-apex[0]) * (vertices[0]-apex[0])
                             + (vertices[1]-apex[1]) * (vertices[1]-apex[1])
                             + (vertices[2]-apex[2]) * (vertices[2]-apex[2]);
        if (this_distance < distance) {
            apex_vertex = i;
            distance = this_distance;
        }
    }
    return apex_vertex;
}

// Convert an array of the form xyzxyz...xyz to a 3D GridFunction
void vertex_vector_to_grid_function(
    mfem::Mesh& mesh,
    double* x,
    mfem::GridFunction *gf)
{
    mfem::FiniteElementCollection *fec = new mfem::H1_FECollection(1, mesh.Dimension());
    mfem::FiniteElementSpace *fespace = new mfem::FiniteElementSpace(&mesh, fec, 3, mfem::Ordering::byVDIM);
    gf->SetSpace(fespace);
    mfem::Vector vals(x, 3*mesh.GetNV());
    gf->SetFromTrueDofs(vals);
}

#ifdef DEBUG
void debug_print_to_file(
    const double *x,
    int n,
    std::string const& dir,
    std::string const& filename)
{
    std::string path = dir + filename;
    std::ofstream fout(path.c_str());
    for (int i = 0; i < n; i++) {
        fout << x[i] << std::endl;
    }
}
#endif
