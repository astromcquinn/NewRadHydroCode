# Making the repository public without Hy Trac's Fortran

*2026-10-26. Hy Trac agreed to NewRadHydroCode being public on the condition
that his original code is removed. This is the plan and the record; the
tree changes are PR "fortran-split", the history step is §3.*

## 1. What is his, and where it went

| Was | Is now |
|---|---|
| `reference/Radhydro/` (the untouched original) | `Radhydro/` in the private repository **`astromcquinn/RadHydro-fortran`** |
| `reference/harness{,_hydro,_pm,_full}/src/` (his modules with the documented patches) and each harness's `Makefile`, run scripts, `README.md`, `PATCHES.md` (which quote his lines), `harness/healpix` (HEALPix 3.11, GPL), `harness/tests` | the same paths in the private repository, so the Makefiles' relative paths (`../Radhydro`, `../harness_common`) still work |
| `reference/harness_common/` (`fft_tools.f90` derived from his, the shared time step copied from `radhydro.f90`, the MKL DFTI interface) | private repository |
| `tests/halo_catalog_ref.f90` (his `read_halo` as a test helper), `tests/fortran_consts/dump_consts.f90` (built against his constants) | `tests/` of the private repository |
| `reference/harness*/outputs*/` (data the harnesses produced) | **stay public**, unchanged paths |

The private repository is the git submodule `reference/fortran` of this
one (`.gitmodules`); it carries the history of those paths (20 commits,
filtered with git-filter-repo) plus a README. With access,
`git submodule update --init reference/fortran` restores the layout the
code and the docs refer to as `reference/fortran/...`.

## 2. What changed in the public tree (PR fortran-split)

* CMake: `RADRAY_FORTRAN_DIR` (default `reference/fortran`);
  `RADRAY_HAVE_FORTRAN_REF` when `harness/Makefile` is there. The live
  gfortran comparisons (`hydro_vs_fortran_live`, `full_vs_fortran_live_*`,
  the halo-catalogue reference program) are defined only then; the
  stored-reference comparisons (`ref_*`, `pm_ref_*`, `hydro_ref_*`,
  `full_ref_*`) are unchanged and still guard the port on every build.
* Workflows: `ci.yml`, `nightly.yml` and `regen_references.yml` populate
  the submodule when the repository secret **`FORTRAN_TOKEN`** is set (a
  fine-grained personal access token with read access to
  `RadHydro-fortran`, "Contents: read"); without it they say so and run
  without the live comparisons. The regeneration workflow runs the
  harnesses from the submodule and stores into the public `outputs*`.
* `CLAUDE.md` rule 1, `README.md`, `reference/README.md`, and the path
  mentions in docs, tools and tests now say `reference/fortran/...`.
  Historical entries in `docs/STATUS.md` and `docs/PLAN.md` keep their old
  paths.
* Verified: the fast suite with the reference absent
  (`-DRADRAY_FORTRAN_DIR=/nonexistent`) and the `fortran`-labelled tests
  with it present, both in the cloud container (§4).

## 3. The history, and the two ways to go public

Every commit since 3a0aeb6 (2026-10-09) holds the Fortran, so a public
repository needs a history without it. `tools/git/scrub_fortran_history.sh`
does the rewrite on a fresh mirror clone (all branches at once) and verifies
that no commit on any ref touches the removed paths and that no blob holds
his module names. Two ways to use it:

**A. A new public repository from the scrubbed history (recommended).**
1. Merge the open thread PRs first, or accept re-creating them (their
   branches are rewritten too, so the PRs can be reopened from the pushed
   branches, but GitHub's PR numbers and discussions start over).
2. Rename this repository to `NewRadHydroCode-private` (GitHub keeps it
   private; its history, issues and PRs stay as the archive).
3. `git clone --no-local --mirror` it, run the scrub, create the new
   **public** `astromcquinn/NewRadHydroCode`, `git push --mirror` into it.
   The old name's redirect stops once the new repository takes the name.
4. Set `FORTRAN_TOKEN`, re-create the open PRs, point every clone at the
   new repository (`git remote set-url` and a fresh `git fetch`; local
   branches must be rebased onto the rewritten commits, or re-cloned).
5. Make nothing else public: `RadHydro-fortran` and the archive stay private.

Nothing of his ever appears in a public history, and nothing has to be
purged from GitHub's caches.

**B. Rewrite in place and force-push, then make this repository public.**
Shorter, keeps the PR numbers, but GitHub keeps unreachable commits and
cached PR views until support purges them, and every clone and open PR
breaks the same way as in A. GitHub's own guidance for removing sensitive
data is this route plus a support ticket; for a repository that has never
been public, A has no residual risk. The user decides.

Either way the rewrite happens once, at a quiet moment: after the current
thread PRs (A, B, C, D of `docs/zoom_turbulence.md` §6) have merged or been
re-based, since their branches are rewritten too.

## 4. Checklist for the user

1. Create the private repository `astromcquinn/RadHydro-fortran` (empty,
   no README); the content and history are ready to push from this session.
2. Create a fine-grained token with "Contents: read" on `RadHydro-fortran`
   and add it to NewRadHydroCode as the Actions secret `FORTRAN_TOKEN`.
3. Merge PR fortran-split (CI runs without the token until step 2).
4. On klone: `git submodule update --init reference/fortran` in the
   checkout (the harness builds and the live tests need it).
5. Decide A or B of §3, and when; then the scrub and the switch to public.
