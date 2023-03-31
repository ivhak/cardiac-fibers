#!/usr/bin/env bash
#
# Copyright (C) 2022 Iver Håkonsen
#
# cardiac-fibers is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with
# cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
#
# Authors: Iver Håkonsen <hakonseniver@yahoo.no

PKGDIR="${HOME}/.local/packages"

MFEM_BUILD_DIR=build-darwin

MFEM_INSTALL_DIR=${PKGDIR}/mfem-4.5

make BUILD_DIR=${MFEM_BUILD_DIR} config \
    MFEM_USE_MPI=YES \
    MPICXX=mpic++ \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L/usr/local/lib -lHYPRE" \
    HYPRE_OPT="-I/usr/local/include" \
    METIS_LIB="-L/usr/local/lib -lmetis" \
    METIS_OPT="-I/usr/local/include"

make BUILD_DIR=${MFEM_BUILD_DIR} -j 4

make BUILD_DIR=${MFEM_BUILD_DIR} check

make BUILD_DIR=${MFEM_BUILD_DIR} install PREFIX=${MFEM_INSTALL_DIR}
