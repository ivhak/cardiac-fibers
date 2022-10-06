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

double timespec_duration(struct timespec t0, struct timespec t1);
void log_timing(std::ostream& out, const char *log_string, double seconds);
std::string basename(std::string const& filename);
std::string remove_extension(std::string const& filename);
void mksubdir(std::string const& subdir);

#endif
