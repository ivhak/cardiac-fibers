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
    out << "[" << std::left << std::setw(12) << log_string << "]: "
        << std::right << std::fixed << std::setw(12) << std::setprecision(6)<< seconds << " s" << std::endl;
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

void mksubdir(std::string const& subdir) {
    const char *cwd = getcwd(NULL, 0);
    std::string path(cwd);
    path += "/";
    path += subdir;
    mkdir(path.c_str(), 0777);
}

