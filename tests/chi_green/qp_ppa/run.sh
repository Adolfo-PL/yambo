#!/bin/bash
#
# License-Identifier: GPL
#
# Copyright (C) 2026 Adolfo Partida
#
# Authors: Adolfo Partida, with AI assistance (Claude Code)
#
# src/qp/QP_ppa_cohsex.F on mock data (stub modules, no Yambo build): the binned retarded
# self-energy (GRetarded) against a copy forced onto Yambo's frequency loop, with RIM-W and
# the BG terminator, 1 and 4 OpenMP threads. Usage: bash tests/chi_green/qp_ppa/run.sh [FC]
set -e
FC=${1:-gfortran}
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../../.." && pwd)
b=$(mktemp -d); trap 'rm -rf "$b"' EXIT
cp "$here"/*.F90 "$here"/*.h "$root/src/modules/mod_Chi_Green.F" "$root/src/qp/QP_ppa_cohsex.F" "$b"/
cd "$b"
sed -e 's/l_PB=QP_retarded_G.and.l_ppa.and..not.have_gpu/l_PB=.FALSE./' \
    -e 's/^subroutine QP_ppa_cohsex(/subroutine QP_ppa_cohsex_direct(/' QP_ppa_cohsex.F > QP_ppa_cohsex_direct.F
grep -q "l_PB=.FALSE." QP_ppa_cohsex_direct.F && grep -q "subroutine QP_ppa_cohsex_direct" QP_ppa_cohsex_direct.F
F="$FC -cpp -fopenmp -ffree-form -ffree-line-length-none -I. -O1 -g -fcheck=bounds"
for f in mocks.F90 mod_Chi_Green.F externals.F90 QP_ppa_cohsex.F QP_ppa_cohsex_direct.F driver.F90; do $F -c $f; done
$F mocks.o mod_Chi_Green.o externals.o QP_ppa_cohsex.o QP_ppa_cohsex_direct.o driver.o -o harness
for t in 1 4; do echo "OMP_NUM_THREADS=$t"; OMP_NUM_THREADS=$t ./harness; done
