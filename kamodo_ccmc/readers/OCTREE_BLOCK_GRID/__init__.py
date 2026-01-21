"""OCTREE_BLOCK_GRID shared library loader (ctypes-based, SpacePy-style).

This module replaces the CFFI-based approach with direct ctypes loading
of the shared library compiled by setup.py.
"""

import ctypes
import os
import sys
import sysconfig
import warnings

# C type aliases for clarity (SpacePy pattern)
c_int_p = ctypes.POINTER(ctypes.c_int)
c_float_p = ctypes.POINTER(ctypes.c_float)


# Define C structures that match the C code
class octree_block(ctypes.Structure):
    """C struct octree_block from octree_block.h"""
    _fields_ = [
        ('refinement_level', ctypes.c_int),
        ('parent_ID', ctypes.c_int),
        ('child_count', ctypes.c_int),
        ('child_IDs', ctypes.c_int * 8),
        ('XMIN', ctypes.c_float),
        ('XMAX', ctypes.c_float),
        ('XCenter', ctypes.c_float),
        ('YMIN', ctypes.c_float),
        ('YMAX', ctypes.c_float),
        ('YCenter', ctypes.c_float),
        ('ZMIN', ctypes.c_float),
        ('ZMAX', ctypes.c_float),
        ('ZCenter', ctypes.c_float),
    ]


def _load_octree_lib():
    """Load OCTREE shared library via ctypes (SpacePy-style)."""
    libdir = os.path.dirname(os.path.abspath(__file__))

    # Try different library naming conventions across platforms
    libnames = {
        'darwin': ['libinterpolate_amrdata.dylib', 'interpolate_amrdata.so'],
        'win32': ['interpolate_amrdata.dll', 'libinterpolate_amrdata.dll'],
        'linux': ['libinterpolate_amrdata.so'],
    }.get(sys.platform, ['libinterpolate_amrdata.so'])

    # Add sysconfig extension suffix variant (SpacePy pattern)
    ext = sysconfig.get_config_var('EXT_SUFFIX')
    if ext is None:
        ext = sysconfig.get_config_var('SO')
    if ext:
        libnames.append('libinterpolate_amrdata' + ext)
        libnames.append('interpolate_amrdata' + ext)

    last_error = None
    for name in libnames:
        libpath = os.path.join(libdir, name)
        if os.path.exists(libpath):
            try:
                lib = ctypes.CDLL(libpath)

                # Set up function signatures (replacing CFFI cdef)
                # int xyz_ranges(int N, float *x, float *y, float *z, ...)
                lib.xyz_ranges.restype = ctypes.c_int
                lib.xyz_ranges.argtypes = [
                    ctypes.c_int,  # N
                    c_float_p,     # x
                    c_float_p,     # y
                    c_float_p,     # z
                    c_float_p,     # XR_BLK
                    c_float_p,     # YR_BLK
                    c_float_p,     # ZR_BLK
                    c_float_p,     # X_BLK
                    c_float_p,     # Y_BLK
                    c_float_p,     # Z_BLK
                    c_float_p,     # box_range
                    c_int_p,       # NX_in
                    c_int_p,       # NY_in
                    c_int_p,       # NZ_in
                    ctypes.c_int,  # positions_in_cell_center
                ]

                # void setup_octree(...)
                lib.setup_octree.restype = None
                lib.setup_octree.argtypes = [
                    ctypes.c_int,                    # N_blks_in
                    c_float_p,                       # xr_blk
                    c_float_p,                       # yr_blk
                    c_float_p,                       # zr_blk
                    ctypes.c_int,                    # MAX_AMRLEVELS
                    c_float_p,                       # box_range
                    ctypes.POINTER(octree_block),   # octree_blocklist_in
                    ctypes.c_int,                    # N_octree_blocks
                    c_int_p,                         # numparents_at_AMRlevel_in
                    c_int_p,                         # block_at_AMRlevel_in
                ]

                # void setup_octree_pointers(...)
                lib.setup_octree_pointers.restype = None
                lib.setup_octree_pointers.argtypes = [
                    ctypes.c_int,                    # MAX_AMRLEVELS_in
                    ctypes.POINTER(octree_block),   # octree_blocklist_in
                    c_int_p,                         # numparents_at_AMRlevel_in
                    c_int_p,                         # block_at_AMRlevel_in
                ]

                # int interpolate_amrdata_multipos(...)
                lib.interpolate_amrdata_multipos.restype = ctypes.c_int
                lib.interpolate_amrdata_multipos.argtypes = [
                    c_float_p,     # xx
                    c_float_p,     # yy
                    c_float_p,     # zz
                    ctypes.c_int,  # npos
                    c_float_p,     # field
                    c_float_p,     # output
                ]

                # float interpolate_amrdata(...)
                lib.interpolate_amrdata.restype = ctypes.c_float
                lib.interpolate_amrdata.argtypes = [
                    ctypes.c_float,  # xx
                    ctypes.c_float,  # yy
                    ctypes.c_float,  # zz
                    c_float_p,       # field
                    ctypes.c_int,    # is_new_position
                ]

                # int trace_fieldline_octree(...)
                lib.trace_fieldline_octree.restype = ctypes.c_int
                lib.trace_fieldline_octree.argtypes = [
                    ctypes.c_float,  # x_start
                    ctypes.c_float,  # y_start
                    ctypes.c_float,  # z_start
                    ctypes.c_float,  # r_end
                    c_float_p,       # bxfield
                    c_float_p,       # byfield
                    c_float_p,       # bzfield
                    c_float_p,       # flx
                    c_float_p,       # fly
                    c_float_p,       # flz
                    c_int_p,         # step_max
                    ctypes.c_float,  # dn
                    ctypes.c_float,  # bdp
                    c_float_p,       # tilt
                    ctypes.c_float,  # spherical_deg2rad
                ]

                return lib

            except OSError as e:
                last_error = e

    # If we get here, we couldn't load the library
    raise RuntimeError(
        f'Cannot load OCTREE library; SWMF-GM reader unavailable. '
        f'Tried: {libnames} in {libdir}. Last error: {last_error}'
    )


# Module-level library loading
try:
    lib = _load_octree_lib()
    OCTREE_AVAILABLE = True
except RuntimeError as e:
    OCTREE_AVAILABLE = False
    lib = None
    warnings.warn(str(e))


# Export the library and availability flag
__all__ = ['lib', 'OCTREE_AVAILABLE', 'octree_block', 'c_int_p', 'c_float_p']
