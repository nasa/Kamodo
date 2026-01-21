"""OpenGGCM shared library loader (ctypes-based, SpacePy-style).

This module loads the OpenGGCM Fortran shared library via ctypes.

COMPILER NOTE: The hidden string length arguments in FORTRAN_FUNCTIONS
assume gfortran's ABI, where string lengths are passed BY VALUE at the
END of the argument list. If compiled with Intel Fortran (ifort), the
ABI is different (lengths immediately follow each string).
"""

import ctypes
import os
import sys
import sysconfig
import warnings
import numpy as np

# Debug mode for array validation (set KAMODO_DEBUG=1 to enable)
KAMODO_DEBUG = os.environ.get('KAMODO_DEBUG', '').lower() in ('1', 'true', 'yes')

# Fortran type aliases (SpacePy pattern for clarity)
int4 = ctypes.c_int32      # Fortran INTEGER (4 bytes)
real4 = ctypes.c_float     # Fortran REAL (4 bytes)
int4_p = ctypes.POINTER(int4)
real4_p = ctypes.POINTER(real4)

# Function signatures dictionary (SpacePy pattern)
# Format: 'name': (restype, [argtypes...])
# Note: Hidden string lengths included for gfortran ABI
FORTRAN_FUNCTIONS = {
    'read_3d_field': (
        None,  # restype (void for subroutines)
        [
            ctypes.c_char_p,   # path3d (character*110)
            real4_p,           # threeDfield
            ctypes.c_char_p,   # var (character*80)
            int4_p,            # fieldsize
            int4_p,            # nx (out)
            int4_p,            # ny (out)
            int4_p,            # nz (out)
            ctypes.c_char_p,   # asciiTime (out, character*80)
            int4,              # hidden: len(path3d)
            int4,              # hidden: len(var)
            int4,              # hidden: len(asciiTime)
        ]
    ),
    'read_2d_field': (
        None,
        [
            ctypes.c_char_p,   # path3d
            real4_p,           # twoDfield
            ctypes.c_char_p,   # var
            int4_p,            # fieldsize
            int4_p,            # nx (out)
            int4_p,            # ny (out)
            ctypes.c_char_p,   # asciiTime (out)
            int4,              # hidden: len(path3d)
            int4,              # hidden: len(var)
            int4,              # hidden: len(asciiTime)
        ]
    ),
    'read_grid_for_vector': (
        None,
        [
            ctypes.c_char_p,   # FileName (character*120)
            ctypes.c_char_p,   # vectorCompName (character*7)
            int4_p,            # nx (out)
            int4_p,            # ny (out)
            int4_p,            # nz (out)
            real4_p,           # gx_v (in/out)
            real4_p,           # gy_v (in/out)
            real4_p,           # gz_v (in/out)
            int4,              # hidden: len(FileName)
            int4,              # hidden: len(vectorCompName)
        ]
    ),
}


def _get_fortran_func(lib, name):
    """Get Fortran function, trying name mangling variations.

    Fortran compilers add suffixes to function names. Common patterns:
    - `read_3d_field_` (gfortran default - single underscore)
    - `read_3d_field__` (some compilers - double underscore)
    - `read_3d_field` (no suffix)
    """
    for suffix in ['_', '__', '']:
        try:
            return getattr(lib, f'{name}{suffix}')
        except AttributeError:
            pass
    raise AttributeError(f'Cannot find Fortran function: {name}')


def _load_openggcm_lib():
    """Load OpenGGCM shared library via ctypes (SpacePy-style)."""
    libdir = os.path.dirname(os.path.abspath(__file__))

    # Try different library naming conventions across platforms
    libnames = {
        'darwin': ['libreadOpenGGCM.dylib', 'libreadOpenGGCM.so'],
        'win32': ['readOpenGGCM.dll', 'libreadOpenGGCM.dll'],
    }.get(sys.platform, ['libreadOpenGGCM.so'])

    # Add sysconfig extension suffix variant (SpacePy pattern)
    ext = sysconfig.get_config_var('EXT_SUFFIX')
    if ext is None:
        ext = sysconfig.get_config_var('SO')
    if ext:
        libnames.append('libreadOpenGGCM' + ext)
        libnames.append('readOpenGGCM' + ext)

    last_error = None
    for name in libnames:
        libpath = os.path.join(libdir, name)
        if os.path.exists(libpath):
            try:
                lib = ctypes.CDLL(libpath)
                _setup_function_signatures(lib)
                return lib
            except OSError as e:
                last_error = e

    raise RuntimeError(
        f'Cannot load OpenGGCM library. Tried: {libnames} in {libdir}. '
        f'Last error: {last_error}'
    )


def _setup_function_signatures(lib):
    """Set up ctypes function signatures for Fortran subroutines.

    Resolves Fortran name mangling once and stores canonical function
    references on the library object. Uses SpacePy's dictionary-based
    pattern for maintainability.

    Note: Hidden string length arguments assume gfortran ABI (lengths
    passed BY VALUE at the END of the argument list). Intel Fortran
    passes lengths immediately after each string argument.

    Raises
    ------
    AttributeError
        If a required function is not found in the library (fail fast).
    """
    for funcname, (restype, argtypes) in FORTRAN_FUNCTIONS.items():
        func = _get_fortran_func(lib, funcname)  # Raises if not found
        func.restype = restype
        func.argtypes = argtypes
        # Store canonical reference (SpacePy pattern)
        setattr(lib, funcname, func)


# ===== High-level wrapper functions =====
# These provide a Python-friendly interface similar to f2py

def read_3d_field(filepath, fielddata, varname):
    """
    Read 3D field from OpenGGCM file.

    Args:
        filepath: Path to the data file (string)
        fielddata: Pre-allocated numpy array, shape (nx, ny, nz), dtype=float32, order='F'
        varname: Variable name to read (string)

    Returns:
        (fielddata, nx, ny, nz, asciitime)
    """
    if lib is None:
        raise RuntimeError('OpenGGCM library not available')

    # Ensure Fortran-contiguous float32
    fielddata = np.asfortranarray(fielddata, dtype=np.float32)

    # Create output variables (passed by reference to Fortran)
    fieldsize = int4(fielddata.size)
    nx = int4(0)
    ny = int4(0)
    nz = int4(0)

    # Fortran character arrays need to be byte strings
    # Fortran expects fixed-length strings, padded with spaces
    path_bytes = filepath.encode('utf-8').ljust(110, b' ')
    var_bytes = varname.encode('utf-8').ljust(80, b' ')
    asciitime_buffer = ctypes.create_string_buffer(80)

    # Use stored reference (SpacePy pattern - no re-resolution)
    func = lib.read_3d_field

    # Call Fortran subroutine
    # Fortran passes strings with hidden length arguments at the end
    func(
        path_bytes,
        fielddata.ctypes.data_as(real4_p),
        var_bytes,
        ctypes.byref(fieldsize),
        ctypes.byref(nx),
        ctypes.byref(ny),
        ctypes.byref(nz),
        asciitime_buffer,
        int4(110),  # length of path3d
        int4(80),   # length of var
        int4(80),   # length of asciiTime
    )

    # Optional debug validation
    if KAMODO_DEBUG:
        assert fielddata.dtype == np.float32, "Array dtype changed after Fortran call!"
        assert fielddata.flags['F_CONTIGUOUS'], "Array lost Fortran contiguity!"

    return fielddata, nx.value, ny.value, nz.value, asciitime_buffer.value


def read_2d_field(filepath, fielddata, varname, nbuffer):
    """
    Read 2D field from OpenGGCM file.

    Args:
        filepath: Path to the data file (string)
        fielddata: Pre-allocated numpy array, 1D, dtype=float32
        varname: Variable name to read (string)
        nbuffer: Size of buffer

    Returns:
        (fielddata, nx, ny, asciitime)
    """
    if lib is None:
        raise RuntimeError('OpenGGCM library not available')

    # Ensure contiguous float32
    fielddata = np.ascontiguousarray(fielddata, dtype=np.float32)

    # Create output variables
    fieldsize = int4(nbuffer)
    nx = int4(0)
    ny = int4(0)

    # Fortran character arrays
    path_bytes = filepath.encode('utf-8').ljust(110, b' ')
    var_bytes = varname.encode('utf-8').ljust(80, b' ')
    asciitime_buffer = ctypes.create_string_buffer(80)

    # Use stored reference (SpacePy pattern - no re-resolution)
    func = lib.read_2d_field

    # Call Fortran subroutine
    func(
        path_bytes,
        fielddata.ctypes.data_as(real4_p),
        var_bytes,
        ctypes.byref(fieldsize),
        ctypes.byref(nx),
        ctypes.byref(ny),
        asciitime_buffer,
        int4(110),  # length of path3d
        int4(80),   # length of var
        int4(80),   # length of asciiTime
    )

    # Optional debug validation
    if KAMODO_DEBUG:
        assert fielddata.dtype == np.float32, "Array dtype changed after Fortran call!"
        assert fielddata.flags['C_CONTIGUOUS'], "Array lost C contiguity!"

    return fielddata, nx.value, ny.value, asciitime_buffer.value


def read_grid_for_vector(gridfile, vecname, gx, gy, gz):
    """
    Read grid coordinates for vector field.

    Args:
        gridfile: Path to grid file (string)
        vecname: Vector component name (string, e.g., 'bx', 'by', 'bz', or ' ' for default)
        gx, gy, gz: Pre-allocated numpy arrays for grid coordinates

    Returns:
        (nx, ny, nz, gx, gy, gz)
    """
    if lib is None:
        raise RuntimeError('OpenGGCM library not available')

    # Ensure contiguous float32
    gx = np.ascontiguousarray(gx, dtype=np.float32)
    gy = np.ascontiguousarray(gy, dtype=np.float32)
    gz = np.ascontiguousarray(gz, dtype=np.float32)

    # Create output dimension variables
    nx = int4(0)
    ny = int4(0)
    nz = int4(0)

    # Fortran character arrays - note vectorCompName is character*7
    file_bytes = gridfile.encode('utf-8').ljust(120, b' ')
    vec_bytes = vecname.encode('utf-8').ljust(7, b' ')

    # Use stored reference (SpacePy pattern - no re-resolution)
    func = lib.read_grid_for_vector

    # Call Fortran subroutine
    func(
        file_bytes,
        vec_bytes,
        ctypes.byref(nx),
        ctypes.byref(ny),
        ctypes.byref(nz),
        gx.ctypes.data_as(real4_p),
        gy.ctypes.data_as(real4_p),
        gz.ctypes.data_as(real4_p),
        int4(120),  # length of FileName
        int4(7),    # length of vectorCompName
    )

    # Optional debug validation
    if KAMODO_DEBUG:
        for arr, name in [(gx, 'gx'), (gy, 'gy'), (gz, 'gz')]:
            assert arr.dtype == np.float32, f"{name} dtype changed after Fortran call!"
            assert arr.flags['C_CONTIGUOUS'], f"{name} lost C contiguity!"

    return nx.value, ny.value, nz.value, gx, gy, gz


# Module-level library loading
try:
    lib = _load_openggcm_lib()
    OPENGGCM_AVAILABLE = True
except RuntimeError as e:
    OPENGGCM_AVAILABLE = False
    lib = None
    warnings.warn(str(e))


# Export
__all__ = [
    'lib', 'OPENGGCM_AVAILABLE',
    'read_3d_field', 'read_2d_field', 'read_grid_for_vector',
    'int4', 'real4', 'int4_p', 'real4_p', 'FORTRAN_FUNCTIONS',
]
