#ifndef LDRB_GPU_HPP
#include <string>
#include "mfem.hpp"

typedef struct {
    int verbose;
    const char *mesh_file;
    std::string mesh_basename;
    const char *output_dir;
    const char *device_config;
    bool paraview;
    mfem::Vector prescribed_apex;
} Options;

typedef enum {
    BASE   = 1,
    EPI    = 2,
    LV_ENDO = 3,
    RV_ENDO = 4,
} MeshAttributes;
#define LDRB_GPU_HPP
#endif
