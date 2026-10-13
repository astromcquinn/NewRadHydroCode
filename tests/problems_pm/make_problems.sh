#!/usr/bin/env bash
# (Re)generate the PM test problems with tools/ics/make_ics.py.
#
#   tests/problems_pm/make_problems.sh [problem ...]
#
# Each problem directory holds params.txt (read by the Fortran harness
# reference/harness_pm and by apps/pm_driver), the initial conditions
# xv_dm.init / u_gas.init (the files io_tools.f90:read_init_conditions reads)
# and ics_info.txt.  The .init files are committed: they are the inputs of
# the stored Fortran reference outputs, so they must not change silently
# (numpy's FFT can differ in the last bit between versions/platforms).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ics="$here/../../tools/ics/make_ics.py"
all="zeldovich_3d growth_3d pancake_x pancake_z_nomf zeldovich_3d_r2 pancake_x_r2 pancake_z_nomf_r2 pancake_x_n32"
problems="${*:-$all}"

params() {  # params <dir> <key=value>...
    local d="$1"; shift
    {
        echo "# PM test problem (tests/problems_pm/make_problems.sh)"
        for kv in "$@"; do echo "${kv%%=*} = ${kv#*=}"; done
    } > "$d/params.txt"
}

for p in $problems; do
    d="$here/$p"
    mkdir -p "$d"
    case "$p" in
    zeldovich_3d)
        # Realistic Zel'dovich ICs: DM and gas follow the same linear field
        # (EH no-wiggle P(k), sigma8 = 0.8, n_s = 0.96), z_init = 300.
        # The code's own time step (calc_timestep: ra_max, amax_dm, cmax_gas).
        python3 "$ics" "$d" --kind zeldovich --Ngas 32 --Ndm 16 --z-init 300 \
            --seed 1 --gas zeldovich --growth full
        params "$d" Ngas=32 Ndm=16 mf=1 z_init=300 nsteps=30 da_sim=0 pk_every=5
        ;;
    growth_3d)
        # Linear growth check: uniform static gas (gravity-only run), DM on
        # the frozen-gas growing mode; P(k) written every step.
        python3 "$ics" "$d" --kind zeldovich --Ngas 32 --Ndm 16 --z-init 300 \
            --seed 2 --gas uniform --growth frozen-gas
        params "$d" Ngas=32 Ndm=16 mf=1 z_init=300 nsteps=30 da_sim=0 pk_every=1
        ;;
    pancake_x)
        # 1D Zel'dovich pancake along x, shell crossing at a = 0.05;
        # code's time step; stops just before crossing.
        python3 "$ics" "$d" --kind pancake --axis x --a-cross 0.05 --Ngas 32 \
            --Ndm 16 --z-init 300 --gas uniform --growth frozen-gas
        params "$d" Ngas=32 Ndm=16 mf=1 z_init=300 nsteps=25 da_sim=0 pk_every=0
        ;;
    pancake_z_nomf)
        # Pancake along z, built WITHOUT -DMF (gravity kicks u_gas instead
        # of v_gas) and with a fixed da_sim = 0.002.
        python3 "$ics" "$d" --kind pancake --axis z --a-cross 0.06 --Ngas 32 \
            --Ndm 16 --z-init 300 --gas uniform --growth frozen-gas
        params "$d" Ngas=32 Ndm=16 mf=0 z_init=300 nsteps=20 da_sim=0.002 pk_every=0
        ;;
    zeldovich_3d_r2)
        # pm_refine = 2 (docs/pm_refine.md): 32^3 PM mesh and 32^3 particles
        # under a 16^3 gas grid; the zeldovich_3d large scales (seed 1).
        python3 "$ics" "$d" --kind zeldovich --Ngas 16 --Ndm 32 --z-init 300 \
            --seed 1 --gas zeldovich --growth full --phases-from-ndm 16
        params "$d" Ngas=16 Ndm=32 pm_refine=2 mf=1 z_init=300 nsteps=30 da_sim=0 pk_every=5
        ;;
    pancake_x_r2)
        # pancake_x on a pm_refine = 2 mesh with Ndm = 2 Ngas
        python3 "$ics" "$d" --kind pancake --axis x --a-cross 0.05 --Ngas 16 \
            --Ndm 32 --z-init 300 --gas uniform --growth frozen-gas
        params "$d" Ngas=16 Ndm=32 pm_refine=2 mf=1 z_init=300 nsteps=25 da_sim=0 pk_every=0
        ;;
    pancake_z_nomf_r2)
        # pancake_z_nomf (gravity kicks u_gas) on a pm_refine = 2 mesh
        python3 "$ics" "$d" --kind pancake --axis z --a-cross 0.06 --Ngas 16 \
            --Ndm 32 --z-init 300 --gas uniform --growth frozen-gas
        params "$d" Ngas=16 Ndm=32 pm_refine=2 mf=0 z_init=300 nsteps=20 da_sim=0.002 pk_every=0
        ;;
    pancake_x_n32)
        # The pm_refine = 1 twin of pancake_x_r2: the same particles
        # (identical xv_dm.init) and the same 32^3 mesh, the gas grid at 32^3.
        # With uniform frozen gas the DM must evolve bit for bit as in
        # pancake_x_r2 (test pm_r2_twin_pancake_x).
        python3 "$ics" "$d" --kind pancake --axis x --a-cross 0.05 --Ngas 32 \
            --Ndm 32 --z-init 300 --gas uniform --growth frozen-gas
        params "$d" Ngas=32 Ndm=32 mf=1 z_init=300 nsteps=25 da_sim=0 pk_every=0
        ;;
    *) echo "unknown problem $p" >&2; exit 2 ;;
    esac
done
