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

#ifndef UTIL_HPP
#define UTIL_HPP
#include <time.h>
#include <fstream>
#include <iomanip>
#include "mfem.hpp"

double timespec_duration(struct timespec t0, struct timespec t1);
void log_timing(std::ostream& out, const char *log_string, double seconds, int ident=0, char override_marker=0);
void log_marker(std::ostream& out, const char *log_string, int ident=0, char override_marker=0);
std::string basename(std::string const& filename);
std::string remove_extension(std::string const& filename);
void mksubdir(std::string const& subdir);

int find_apex_vertex(mfem::Mesh& mesh, mfem::Vector& apex);
void vertex_vector_to_grid_function(mfem::Mesh& mesh, double* x, mfem::GridFunction *gf);
void save_solution(mfem::GridFunction *x, std::string const& dir, std::string const& base_name, std::string const& suffix);
#ifdef MFEM_USE_MPI
void save_solution(mfem::ParGridFunction *x, std::string const& dir, std::string const& base_name, std::string const& suffix, int rank);
#endif
#ifdef DEBUG
void debug_print_to_file(const double* x, int n, std::string const& dir, std::string const& filename);
#endif

#endif
