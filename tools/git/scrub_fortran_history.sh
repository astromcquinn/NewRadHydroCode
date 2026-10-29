#!/bin/bash
# Rewrite the history of NewRadHydroCode without Hy Trac's Fortran (the step
# before the repository is made public; docs/fortran_split.md). Works on a
# FRESH clone made with --no-local (git filter-repo refuses anything else),
# removes every path that ever held the Fortran or quoted it, and then
# verifies that no commit on any branch still touches them. The result is
# pushed by hand to the new public repository (never force-pushed over the
# private one without the user's say-so: every clone, branch and open PR
# would have to be re-based).
#
#   git clone --no-local --mirror https://github.com/astromcquinn/NewRadHydroCode.git nrhc-scrub
#   tools/git/scrub_fortran_history.sh nrhc-scrub
#
# A mirror clone carries every branch and tag, so the rewrite applies to all
# of them at once (the threads' branches included). pip install git-filter-repo.
set -eu
repo="${1:?usage: scrub_fortran_history.sh <fresh mirror clone>}"
cd "$repo"
git filter-repo --force --invert-paths \
  --path Radhydro --path harness --path harness_common --path harness_full --path harness_hydro --path harness_pm \
  --path reference/Radhydro \
  --path reference/harness/src \
  --path reference/harness/tests \
  --path reference/harness/PATCHES.md \
  --path reference/harness_hydro/src \
  --path reference/harness_hydro/PATCHES.md \
  --path reference/harness_pm/src \
  --path reference/harness_pm/PATCHES.md \
  --path reference/harness_full/src \
  --path reference/harness_full/PATCHES.md \
  --path reference/harness_common \
  --path tests/halo_catalog_ref.f90 \
  --path tests/fortran_consts/dump_consts.f90
# the Makefiles, run scripts and HEALPix (GPL) stayed in the history: ours or
# redistributable; the PATCHES.md files quote his lines and go. The root-level
# Radhydro/ is where the original sat in the first commits, before reference/
# (found by the verification below on 2026-10-26: the first pass missed it).
echo "== verification: commits on any ref that still touch the removed paths (must be none)"
n=$(git log --all --format=%h -- Radhydro reference/Radhydro 'reference/harness*/src' reference/harness_common \
      'reference/harness*/PATCHES.md' tests/halo_catalog_ref.f90 tests/fortran_consts/dump_consts.f90 | wc -l)
echo "commits: $n"
[ "$n" = 0 ] || { echo "SCRUB INCOMPLETE"; exit 1; }
echo "== Fortran sources left anywhere in the history (expected: only HEALPix under reference/harness/healpix):"
git log --all --name-only --format= -- '*.f90' '*.F90' '*.f' | sort -u | sed 's#/[^/]*$##' | sort -u
echo "== a grep of every blob for Hy Trac's module names (must be empty):"
git rev-list --all --objects | awk '{print $1}' | git cat-file --batch-check='%(objecttype) %(objectname)' | awk '$1=="blob"{print $2}' \
  | while read b; do git cat-file -p "$b" | grep -l -m1 -E "^ *module (ray_tracing_tools|rt_tools|hydro_tools|pm_tools|halo_tools|radhydro_global)\b" >/dev/null 2>&1 && echo "blob $b"; done || true
git reflog expire --expire=now --all && git gc --prune=now --aggressive -q
echo "== done: $(git rev-list --all --count) commits, $(du -sh .git 2>/dev/null | cut -f1)"
