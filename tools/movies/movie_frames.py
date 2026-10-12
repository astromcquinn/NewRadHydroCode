"""Readers for apps/radhydro movie frames: full cubes and slabs (docs/movies.md).

A frame field is either

  movie/frame_NNNNN_<field>.bin    float32 (N,N,N), Fortran order (x fastest):
                                   the default (movie_slab = 0); or
  movie/frame_NNNNN_<field>.slab   movie_slab = t > 0: for each axis in
                                   `axes` of movie/slab.txt in turn, the block
                                   of cells start .. start+t-1 (0-based) along
                                   that axis and all N along the others,
                                   float32 in Fortran order (ncomp, nx, ny, nz)
                                   with n_axis = t; ncomp = 1 except gamma (3)
                                   and ngam (Nfreq).

movie/slab.txt (key = value; '#' comments) holds N, thick, axes and start.
The whole-box means x_HII,V and x_HII,M are in movie/frames.txt either way.
A run restarted with movie_slab changed from 0 may hold both kinds; each
frame file is read by its extension.
"""
import os

import numpy as np

AXES = "xyz"
_INFO = {}


def slab_info(movie_dir):
    """dict(N, thick, axes, start) from movie_dir/slab.txt, or None."""
    movie_dir = os.path.abspath(movie_dir)
    if movie_dir not in _INFO:
        p = os.path.join(movie_dir, "slab.txt")
        info = None
        if os.path.exists(p):
            d = {}
            for line in open(p):
                line = line.split("#")[0]
                if "=" in line:
                    k, v = line.split("=", 1)
                    d[k.strip()] = v.strip()
            info = dict(N=int(d["N"]), thick=int(d["thick"]), axes=d["axes"],
                        start=int(d["start"].split()[0]))
        _INFO[movie_dir] = info
    return _INFO[movie_dir]


def run_slab_info(rundir):
    return slab_info(os.path.join(rundir, "movie"))


def frame_file(rundir, frame, field):
    """Path of a frame's field (.bin preferred, else .slab), or None."""
    for ext in (".bin", ".slab"):
        p = os.path.join(rundir, "movie", "frame_%05d_%s%s" % (frame, field, ext))
        if os.path.exists(p):
            return p
    return None


def is_slab(path):
    return path.endswith(".slab")


def info_of(path):
    info = slab_info(os.path.dirname(path))
    if info is None:
        raise ValueError("%s: slab frame without movie/slab.txt" % path)
    return info


def cube_N(path):
    """N of a float32 N^3 cube file (None if the size is not a cube)."""
    nb = os.path.getsize(path)
    n = round((nb // 4) ** (1 / 3))
    return n if n ** 3 * 4 == nb else None


def grid_N(path):
    """N of the grid behind a frame file (cube or slab)."""
    return info_of(path)["N"] if is_slab(path) else cube_N(path)


def slab_block(path, axis):
    """Memory-mapped block of a .slab file normal to axis (0, 1, 2): array
    of shape (N, N, N) with n_axis = thick, Fortran order, holding cells
    start .. start+thick-1 along axis."""
    info = info_of(path)
    N, t, axes = info["N"], info["thick"], info["axes"]
    if AXES[axis] not in axes:
        raise ValueError("%s holds slabs normal to %s only, not %s" % (path, ",".join(axes), AXES[axis]))
    nblk = N * N * t
    ncomp = os.path.getsize(path) // (4 * nblk * len(axes))
    if ncomp != 1 or os.path.getsize(path) != 4 * nblk * len(axes):
        raise ValueError("%s: %d B is not %d slab(s) of one float32 component"
                         % (path, os.path.getsize(path), len(axes)))
    shape = [N, N, N]
    shape[axis] = t
    k = axes.index(AXES[axis])
    return np.memmap(path, dtype="<f4", mode="r", offset=4 * nblk * k, shape=tuple(shape), order="F")


def plane_error(path, axis, i0, thick):
    """None if cells i0 .. i0+thick-1 along axis are in the file, else why not."""
    if not is_slab(path):
        return None
    info = info_of(path)
    s0, t = info["start"], info["thick"]
    if AXES[axis] not in info["axes"]:
        return ("the slab frames hold slabs normal to %s only (movie_slab_axes); use --axis %s"
                % (",".join(info["axes"]), info["axes"][0]))
    if i0 < s0 or i0 + thick > s0 + t:
        return ("the slab frames hold cells %d..%d along %s (movie_slab = %d); cells %d..%d "
                "requested: choose --pos/--thick inside the slab" % (s0, s0 + t - 1, AXES[axis], t,
                                                                     i0, i0 + thick - 1))
    return None


def default_i0(N, thick, info):
    """First cell of the default plane: N//2 for full cubes (as before);
    with slab frames the plane of thick cells centred on N/2, kept inside
    the slab (thick = 1: also N//2)."""
    if info is None:
        return N // 2
    i0 = N // 2 if thick == 1 else N // 2 - thick // 2
    return int(min(max(i0, info["start"]), info["start"] + info["thick"] - thick))


def _cells(a, axis, i0, thick):
    """np.take(a, [(i0 + k) % n for k < thick], axis) of a memory-mapped
    array, as the same C-contiguous float32 array, but read by basic
    slicing: np.take first copies the whole (Fortran-ordered) file to a
    C-ordered array, i.e. reads every cell (0.1 s at 256^3 from the page
    cache, minutes for 512^3 frames on GPFS with many processes), while
    a slice reads only the pages of those cells."""
    n = a.shape[axis]
    parts, i, left = [], i0 % n, thick
    while left > 0:
        k = min(left, n - i)
        sl = [slice(None)] * 3
        sl[axis] = slice(i, i + k)
        parts.append(np.ascontiguousarray(np.asarray(a[tuple(sl)])))
        left -= k
        i = 0
    return parts[0] if len(parts) == 1 else np.concatenate(parts, axis=axis)


def read_cells(path, axis, i0, thick, N=None):
    """Cells i0 .. i0+thick-1 (periodic for cubes) along axis of a frame
    file (cube or slab), float32, C-contiguous, shape (N, N, N) with
    n_axis = thick: equal to np.take(cube, ..., axis) of the full cube."""
    if is_slab(path):
        err = plane_error(path, axis, i0, thick)
        if err:
            raise ValueError("%s: %s" % (path, err))
        return _cells(slab_block(path, axis), axis, i0 - info_of(path)["start"], thick)
    N = N or cube_N(path)
    a = np.memmap(path, dtype="<f4", mode="r", shape=(N, N, N), order="F")
    return _cells(a, axis, i0, thick)


def read_plane(path, axis, i0, thick, N=None):
    """Mean over cells i0 .. i0+thick-1 (periodic for cubes) along axis;
    2-D float64 array indexed [first remaining axis, second remaining
    axis]. Reads only those cells."""
    return read_cells(path, axis, i0, thick, N).astype(np.float64).mean(axis=axis)


def axis_profile(path, axis, N=None):
    """Mean over the planes perpendicular to axis, as a function of the
    cell along axis (length N). For a slab file the mean is over the first
    stored slab that contains the axis direction (not the whole box); the
    second value is a note saying so ('' for full cubes)."""
    if is_slab(path):
        info = info_of(path)
        for a in info["axes"]:
            b = AXES.index(a)
            if b != axis:
                blk = np.asarray(slab_block(path, b), dtype=np.float64)
                other = tuple(i for i in range(3) if i != axis)
                return blk.mean(axis=other), "slab normal to %s only" % a
        raise ValueError("%s: no stored slab contains the %s direction" % (path, AXES[axis]))
    N = N or cube_N(path)
    a = np.memmap(path, dtype="<f4", mode="r", shape=(N, N, N), order="F")
    other = tuple(i for i in range(3) if i != axis)
    return np.asarray(a, dtype=np.float64).mean(axis=other), ""
