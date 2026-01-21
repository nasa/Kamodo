"""Tri2D shared library loader (ctypes-based, SpacePy-style).

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


def _load_tri2d_lib():
    """Load Tri2D shared library via ctypes (SpacePy-style)."""
    libdir = os.path.dirname(os.path.abspath(__file__))

    # Try different library naming conventions across platforms
    libnames = {
        'darwin': ['libinterpolate_tri2d.dylib', 'interpolate_tri2d.so'],
        'win32': ['interpolate_tri2d.dll', 'libinterpolate_tri2d.dll'],
        'linux': ['libinterpolate_tri2d.so'],
    }.get(sys.platform, ['libinterpolate_tri2d.so'])

    # Add sysconfig extension suffix variant (SpacePy pattern)
    ext = sysconfig.get_config_var('EXT_SUFFIX')
    if ext is None:
        ext = sysconfig.get_config_var('SO')
    if ext:
        libnames.append('libinterpolate_tri2d' + ext)
        libnames.append('interpolate_tri2d' + ext)

    last_error = None
    for name in libnames:
        libpath = os.path.join(libdir, name)
        if os.path.exists(libpath):
            try:
                lib = ctypes.CDLL(libpath)

                # Set up function signatures (replacing CFFI cdef)
                # void setup_tri_pointers(...)
                lib.setup_tri_pointers.restype = None
                lib.setup_tri_pointers.argtypes = [
                    c_float_p,  # tri_x_in
                    c_float_p,  # tri_y_in
                    c_int_p,    # tri_nx_in
                    c_int_p,    # tri_ny_in
                    c_int_p,    # n_tri_in
                    c_int_p,    # tri_vertices_in
                    c_int_p,    # tri_ij_in
                    c_int_p,    # tri_neighbors_in
                    c_float_p,  # tri_xmin_sg_in
                    c_float_p,  # tri_ymin_sg_in
                    c_float_p,  # tri_dx_sg_in
                    c_float_p,  # tri_dy_sg_in
                    c_int_p,    # tri_nx_sg_in
                    c_int_p,    # tri_ny_sg_in
                    c_int_p,    # tri_start_index_sg_in
                    c_float_p,  # tri_z_in
                    c_int_p,    # tri_nz_in
                ]

                # int interpolate_tri2d_plus_1d_multipos(...)
                lib.interpolate_tri2d_plus_1d_multipos.restype = ctypes.c_int
                lib.interpolate_tri2d_plus_1d_multipos.argtypes = [
                    c_float_p,     # xx
                    c_float_p,     # yy
                    c_float_p,     # zz
                    ctypes.c_int,  # npos
                    c_float_p,     # field
                    c_float_p,     # output
                    ctypes.c_int,  # data_in_center
                ]

                return lib

            except OSError as e:
                last_error = e

    # If we get here, we couldn't load the library
    raise RuntimeError(
        f'Cannot load Tri2D library; GAMERA-GM reader unavailable. '
        f'Tried: {libnames} in {libdir}. Last error: {last_error}'
    )


# Module-level library loading
try:
    lib = _load_tri2d_lib()
    TRI2D_AVAILABLE = True
except RuntimeError as e:
    TRI2D_AVAILABLE = False
    lib = None
    warnings.warn(str(e))


# Export the library and availability flag
__all__ = ['lib', 'TRI2D_AVAILABLE', 'c_int_p', 'c_float_p']
