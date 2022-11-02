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

#ifdef HIP_TRACE
#include <roctx.h>
#endif

#include "util.hpp"

void util::tracing::roctx_range_push(const char *s) {
#ifdef HIP_TRACE
    roctxRangePush(s);
#endif
}

void util::tracing::roctx_range_pop(void) {
#ifdef HIP_TRACE
    roctxRangePop();
#endif
}

static const char level_markers[] = {'-', '*', '+'};
#define INDENT_WIDTH 2
#define LOG_LEFT_COL_WIDTH 30

void util::logging::timestamp(
    std::ostream& out,
    std::string const& log_string,
    double seconds,
    const int indent,
    char override_marker)
{
    assert(indent >= 0 && indent <= 2);
    char level_marker = level_markers[indent];

    if (override_marker)
        level_marker = override_marker;

    std::string indent_chars = std::string(INDENT_WIDTH*indent, ' ');
    out << indent_chars
        << level_marker
        << " "
        << std::left << std::setw(LOG_LEFT_COL_WIDTH-INDENT_WIDTH*indent) << log_string
        << std::right << std::fixed << std::setw(12)
        << std::setprecision(6)<< seconds << " s" << std::endl;
}

void util::logging::marker(
    std::ostream& out,
    std::string const& log_string,
    const int indent,
    char override_marker)
{
    assert(indent >= 0 && indent <= 2);
    char level_marker = level_markers[indent];

    std::string indent_chars = std::string(INDENT_WIDTH*indent, ' ');

    if (override_marker)
        level_marker = override_marker;

    out << indent_chars
        << level_marker
        << " "
        << log_string
        << std::endl;
}

void util::logging::info(
    std::ostream& out,
    std::string const& log_string)
{
    out << "[INFO]: " << log_string << std::endl;
}

double util::timing::duration(
    struct timespec t0,
    struct timespec t1)
{
    return (t1.tv_sec - t0.tv_sec) +
        (t1.tv_nsec - t0.tv_nsec) * 1e-9;
}

void util::timing::tick(struct timespec *t)
{
    clock_gettime(CLOCK_MONOTONIC, t);
}



std::string util::fs::basename(std::string const& path)
{
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string util::fs::remove_extension(std::string const& filename)
{
    size_t index_of_extension = filename.find_last_of(".");
    return filename.substr(0, index_of_extension);
}

void util::fs::mksubdir(std::string const& subpath)
{
    const char *cwd = getcwd(NULL, 0);
    std::string cwdpath(cwd);
    std::string path;

    size_t pos = 0;
    std::string subdir;
    while ((pos = subpath.find("/", pos)) != std::string::npos) {
        subdir = subpath.substr(0, pos);
        path = cwdpath + "/" + subdir;
        mkdir(path.c_str(), 0777);
        pos += subdir.length();
    }
    path = cwdpath + "/" + subpath;
    mkdir(path.c_str(), 0777);
}

// Save a solution (in form of a GridFunction) to a file named "<prefix><suffix>".
void util::save::save_solution(
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

#ifdef MFEM_USE_MPI
// Save a solution (in form of a GridFunction) to a file named "<prefix><suffix>".
void util::save::save_solution(
    mfem::ParGridFunction *x,
    std::string const& dir,
    std::string const& base_name,
    std::string const& suffix,
    int rank)
{
    std::ostringstream filename;
    filename << dir << "/" << base_name << suffix << "."
             << std::setfill('0') << std::setw(6) << rank;
    std::ofstream x_ofs(filename.str().c_str());
    x_ofs.precision(8);
    x->Save(x_ofs);
}

void util::save::save_mesh(
    mfem::ParMesh *pmesh,
    std::string const& dir,
    std::string const& base_name,
    int rank)
{
    std::ostringstream mesh_out;
    mesh_out << dir << "/" << base_name << ".mesh."
             << std::setfill('0') << std::setw(6) << rank;
    std::ofstream mesh_ofs(mesh_out.str().c_str());
    mesh_ofs.precision(8);
    pmesh->Print(mesh_ofs);
}
#endif

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
