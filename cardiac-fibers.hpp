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

#ifndef CARDIAC_FIBERS_HPP
#include <string>
#include "mfem.hpp"

typedef struct {
    int verbose;
    const char *mesh_file;
    const char *output_dir;
    const char *device_config;
    bool time_to_file;
    bool save_paraview;
    bool save_mfem;
    mfem::Vector prescribed_apex;
    double itol;
    double alpha_endo;
    double alpha_epi;
    double beta_endo;
    double beta_epi;
    int base_id;
    int epi_id;
    int lv_id;
    int rv_id;
    int uniform_refinement;
#if defined(MFEM_USE_HIP) || defined(MFEM_USE_CUDA)
    bool gpu_tuned_amg;
#endif
    bool use_dg;
} Options;

#define LDRB_GPU_HPP
#endif
