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

## 3b. Done 2026-10-26: route A, the public repository

* PR #61 (the tree) merged in the private repository on the local
  verification (Actions were blocked by the budget).
* The user renamed the repository to **`NewRadHydroCode-private`** (the
  archive: full history, issues, PRs 1–61) and created the new, public
  **`astromcquinn/NewRadHydroCode`**.
* A mirror clone of the private repository was scrubbed with
  `tools/git/scrub_fortran_history.sh` (two passes: the script's paths,
  then the root-level `Radhydro/` of the first commits, now in the
  script), verified (no commit on any ref touches the removed paths; the
  only Fortran left is HEALPix; no blob holds his module names) and pushed
  with `--mirror`: 42 branches, 744 commits, `main` at db2f4b4 (the
  rewritten merge of #61). Pull-request refs were dropped before the push.
* The open PRs were re-opened on the new repository with their
  descriptions: private #60 (`zoom-ics`) is public **#1**, private #58
  (`rt-suite-512`) is public **#2**. PR numbers start over there.
* **Every clone and worktree made before the move must be re-pointed**;
  the branches have the same content under new commit ids:

  ```
  git remote set-url origin https://github.com/astromcquinn/NewRadHydroCode.git
  git fetch origin
  git rebase --onto origin/<branch> <old base> <branch>   # or, with no local commits: git reset --hard origin/<branch>
  git submodule update --init reference/fortran           # with access to RadHydro-fortran
  ```

  A push of a pre-move branch is refused by the new repository's history
  (non-fast-forward) rather than silently reintroducing the Fortran; the
  private repository still accepts them, so keep its remote out of
  everyday work.
* Actions on the public repository are free on the standard runners, so
  the budget block of 2026-10-25 no longer applies there.

## 4. Checklist for the user

1. ~~Create the private repository `astromcquinn/RadHydro-fortran`~~ done
   (`main` at 6e7e31b).
2. Create a fine-grained token with "Contents: read" on `RadHydro-fortran`
   and add it to the public NewRadHydroCode as the Actions secret
   `FORTRAN_TOKEN` (until then CI skips the live comparisons).
3. ~~Merge PR fortran-split~~ done (#61 in the private repository).
4. On klone: re-point the checkout (§3b) and
   `git submodule update --init reference/fortran`.
5. ~~Decide A or B~~ A, done (§3b). Tell the running threads to re-point
   (§3b); their PRs are #1 and #2 on the new repository.
