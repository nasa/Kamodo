"""
Custom build commands for kamodo-ccmc

This module provides custom setuptools commands to automatically compile
C and Fortran extensions during pip install, following SpacePy's proven approach.

Extensions compiled:
- OCTREE_BLOCK_GRID (C + CFFI) - for SWMF-GM model reader
- Tri2D (C + CFFI) - for GAMER-AM model reader
- OpenGGCM (Fortran + f2py) - for OpenGGCM model reader
"""

import os
import sys
import subprocess
import warnings
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install


def compile_extensions(skip_extensions=False):
    """
    Compile all C and Fortran extensions.

    This is the main entry point called by custom install/develop commands.
    """
    if skip_extensions:
        warnings.warn("Skipping C/Fortran extension compilation (--skip-extensions specified)")
        return

    # Check if KAMODO_RELEASE is set (strict mode - fail on errors)
    release_build = bool(os.environ.get("KAMODO_RELEASE", False))

    print("\n" + "="*70)
    print("Kamodo: Compiling C and Fortran extensions")
    print("="*70 + "\n")

    # Track which extensions succeeded
    octree_ok = _compile_octree_extension(release_build)
    tri2d_ok = _compile_tri2d_extension(release_build)
    openggcm_ok = _compile_openggcm_fortran(release_build)

    # Summary
    print("\n" + "="*70)
    print("Extension compilation summary:")
    print(f"  OCTREE_BLOCK_GRID (C):  {'OK' if octree_ok else 'FAILED'}")
    print(f"  Tri2D (C):              {'OK' if tri2d_ok else 'FAILED'}")
    print(f"  OpenGGCM (Fortran):     {'OK' if openggcm_ok else 'FAILED'}")

    if not (octree_ok and tri2d_ok and openggcm_ok):
        msg = (
            "\nWARNING: Some extensions failed to compile. kamodo-ccmc will still "
            "install, but some model readers will be unavailable:\n"
        )
        if not octree_ok:
            msg += "  - SWMF-GM reader (requires OCTREE_BLOCK_GRID)\n"
        if not tri2d_ok:
            msg += "  - GAMER-AM reader (requires Tri2D)\n"
        if not openggcm_ok:
            msg += "  - OpenGGCM reader (requires OpenGGCM Fortran)\n"
        msg += "\nTo fix: Install gcc/gfortran and reinstall kamodo-ccmc\n"
        warnings.warn(msg)

    print("="*70 + "\n")


class install(_install):
    """Custom install command that compiles extensions after installation"""

    user_options = _install.user_options + [
        ('skip-extensions', None, 'Skip compiling C/Fortran extensions'),
    ]

    def initialize_options(self):
        super().initialize_options()
        self.skip_extensions = False

    def run(self):
        # Run the standard install first
        super().run()
        # Then compile extensions
        compile_extensions(self.skip_extensions)


class develop(_develop):
    """Custom develop command that compiles extensions after installation"""

    user_options = _develop.user_options + [
        ('skip-extensions', None, 'Skip compiling C/Fortran extensions'),
    ]

    def initialize_options(self):
        super().initialize_options()
        self.skip_extensions = False

    def run(self):
        # Run the standard develop first
        super().run()
        # Then compile extensions
        compile_extensions(self.skip_extensions)


class build_ext(_build_ext):
    """Custom build_ext - kept for backwards compatibility and explicit builds"""

    user_options = _build_ext.user_options + [
        ('skip-extensions', None, 'Skip compiling C/Fortran extensions'),
    ]

    def initialize_options(self):
        super().initialize_options()
        self.skip_extensions = False

    def run(self):
        """Run standard build_ext, then compile our extensions"""
        super().run()
        if not self.skip_extensions:
            compile_extensions(skip_extensions=False)


def _compile_octree_extension(release_build=False):
    """Compile OCTREE_BLOCK_GRID C extension using CFFI"""
    print("\n--- Compiling OCTREE_BLOCK_GRID (C + CFFI) ---")

    # Path to the OCTREE_BLOCK_GRID directory
    octree_dir = os.path.join(
        os.path.dirname(__file__),
        'readers', 'OCTREE_BLOCK_GRID'
    )

    if not os.path.exists(octree_dir):
        msg = f"OCTREE_BLOCK_GRID directory not found: {octree_dir}"
        if release_build:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False

    # Save current directory
    orig_dir = os.getcwd()

    try:
        os.chdir(octree_dir)

        # Import and run the CFFI builder
        try:
            # Try importing the builder module
            sys.path.insert(0, octree_dir)
            from interpolate_amrdata_extension_build import build_extension

            # Run the build
            success = build_extension()

            if success:
                print("OCTREE_BLOCK_GRID compilation successful")
                return True
            else:
                msg = "OCTREE_BLOCK_GRID compilation failed"
                if release_build:
                    raise RuntimeError(msg)
                warnings.warn(msg)
                return False

        except ImportError as e:
            # Fallback: Try running as script (old way)
            print(f"Could not import builder module ({e}), trying script execution...")
            retval = subprocess.call(
                [sys.executable, 'interpolate_amrdata_extension_build.py'],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )

            if retval == 0:
                print("OCTREE_BLOCK_GRID compilation successful (script mode)")
                return True
            else:
                msg = f"OCTREE_BLOCK_GRID compilation failed with exit code {retval}"
                if release_build:
                    raise RuntimeError(msg)
                warnings.warn(msg)
                return False

    except Exception as e:
        msg = f"OCTREE_BLOCK_GRID compilation error: {e}"
        if release_build:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False
    finally:
        os.chdir(orig_dir)
        if octree_dir in sys.path:
            sys.path.remove(octree_dir)


def _compile_tri2d_extension(release_build=False):
    """Compile Tri2D C extension using CFFI"""
    print("\n--- Compiling Tri2D (C + CFFI) ---")

    # Path to the Tri2D directory
    tri2d_dir = os.path.join(
        os.path.dirname(__file__),
        'readers', 'Tri2D'
    )

    if not os.path.exists(tri2d_dir):
        msg = f"Tri2D directory not found: {tri2d_dir}"
        if release_build:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False

    # Save current directory
    orig_dir = os.getcwd()

    try:
        os.chdir(tri2d_dir)

        # Import and run the CFFI builder
        try:
            # Try importing the builder module
            sys.path.insert(0, tri2d_dir)
            from interpolate_tri2d_extension_build import build_extension

            # Run the build
            success = build_extension()

            if success:
                print("Tri2D compilation successful")
                return True
            else:
                msg = "Tri2D compilation failed"
                if release_build:
                    raise RuntimeError(msg)
                warnings.warn(msg)
                return False

        except ImportError as e:
            # Fallback: Try running as script (old way)
            print(f"Could not import builder module ({e}), trying script execution...")
            retval = subprocess.call(
                [sys.executable, 'interpolate_tri2d_extension_build.py'],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )

            if retval == 0:
                print("Tri2D compilation successful (script mode)")
                return True
            else:
                msg = f"Tri2D compilation failed with exit code {retval}"
                if release_build:
                    raise RuntimeError(msg)
                warnings.warn(msg)
                return False

    except Exception as e:
        msg = f"Tri2D compilation error: {e}"
        if release_build:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False
    finally:
        os.chdir(orig_dir)
        if tri2d_dir in sys.path:
            sys.path.remove(tri2d_dir)


def _compile_openggcm_fortran(release_build=False):
    """Compile OpenGGCM Fortran extension using f2py"""
    print("\n--- Compiling OpenGGCM (Fortran + f2py) ---")

    # Path to the OpenGGCM directory
    openggcm_dir = os.path.join(
        os.path.dirname(__file__),
        'readers', 'OpenGGCM'
    )

    if not os.path.exists(openggcm_dir):
        msg = f"OpenGGCM directory not found: {openggcm_dir}"
        if release_build:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False

    # Save current directory
    orig_dir = os.getcwd()

    try:
        os.chdir(openggcm_dir)

        # Check if f2py is available
        try:
            # Try numpy.f2py first (modern way)
            import numpy.f2py
            f2py_cmd = [sys.executable, '-m', 'numpy.f2py']
        except ImportError:
            # Fall back to f2py command (older numpy versions)
            f2py_cmd = ['f2py']

        # Build command: f2py -c -m readOpenGGCM read_b_grids.f readmagfile3d.f
        cmd = f2py_cmd + [
            '-c', '-m', 'readOpenGGCM',
            'read_b_grids.f', 'readmagfile3d.f'
        ]

        print(f"Running: {' '.join(cmd)}")

        # Run f2py
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )

        if result.returncode == 0:
            print("OpenGGCM compilation successful")
            return True
        else:
            msg = (
                f"OpenGGCM compilation failed with exit code {result.returncode}\n"
                f"Output: {result.stdout[:500] if result.stdout else 'None'}\n"
                f"Hint: Install gfortran and meson, and ensure numpy is installed"
            )
            if release_build:
                raise RuntimeError(msg)
            warnings.warn(msg)
            return False

    except Exception as e:
        msg = f"OpenGGCM compilation error: {e}"
        if release_build:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False
    finally:
        os.chdir(orig_dir)
