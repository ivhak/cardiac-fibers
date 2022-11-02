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

#ifndef LDRB_GPU_HPP
#include <string>
#include "mfem.hpp"

typedef struct {
    int verbose;
    const char *mesh_file;
    std::string mesh_basename;
    const char *output_dir;
    const char *device_config;
    const char *time_file;
    bool paraview;
    mfem::Vector prescribed_apex;
    double alpha_endo;
    double alpha_epi;
    double beta_endo;
    double beta_epi;
    int base_id;
    int epi_id;
    int lv_id;
    int rv_id;
    bool geom_has_rv;
    int solver;
    bool gpu_tuned_amg;
} Options;

#define LDRB_GPU_HPP
#endif
