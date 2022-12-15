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

#ifndef UTIL_HPP
#define UTIL_HPP
#include <time.h>
#include <fstream>
#include <iomanip>
#include "mfem.hpp"

#define INDENT_WIDTH 2
#define LOG_LEFT_COL_WIDTH 48

namespace util {
    namespace tracing {
        // Wrappers around roctxRangePush/Pop that are no-ops if tracing is not enabled
        void roctx_range_push(const char *s);
        void roctx_range_pop();
    }

    namespace logging {
        void timestamp(std::ostream& out, std::string const& log_string, double seconds,
                       int ident=0, char override_marker=0, int colwidth=LOG_LEFT_COL_WIDTH);
        void marker(std::ostream& out, std::string const& log_string,
                int ident=0, char override_marker=0);
        void info(std::ostream& out, std::string const& log_string);
    }

    namespace timing {
        void tick(struct timespec *t);
        double duration(struct timespec t0, struct timespec t1);
    }

    namespace fs {
        std::string basename(std::string const& filename);
        std::string remove_extension(std::string const& filename);
        void mksubdir(std::string const& subdir);
    }

    namespace save {
        void save_solution(mfem::GridFunction *x, std::string const& dir,
                           std::string const& base_name, std::string const& suffix);
#ifdef MFEM_USE_MPI
        void save_solution(mfem::ParGridFunction *x, std::string const& dir,
                           std::string const& base_name, std::string const& suffix,
                           int rank);
        void save_mesh(mfem::ParMesh *pmesh, std::string const& dir,
                       std::string const& base_name, int rank);
#endif
    }
}

#ifdef DEBUG
void debug_print_to_file(const double* x, int n, std::string const& dir, std::string const& filename);
#endif

#endif
