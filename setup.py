from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import subprocess
import sys
import os
import shutil
import glob
import warnings

# Use CFFI's native setuptools integration to compile C extensions during wheel build.
# This tells setuptools to invoke each ffibuilder.compile() during the build phase.
# The format is "path/to/build_script.py:ffibuilder_object_name"
cffi_modules = [
    "kamodo_ccmc/readers/OCTREE_BLOCK_GRID/interpolate_amrdata_extension_build.py:ffibuilder",
    "kamodo_ccmc/readers/Tri2D/interpolate_tri2d_extension_build.py:ffibuilder",
]


class build_ext(_build_ext):
    """Custom build_ext that compiles OpenGGCM Fortran extension via f2py.

    This follows SpacePy's approach: register an Extension to trigger build_ext,
    then do the actual compilation in run() via subprocess.
    """

    def run(self):
        # First, run the standard build_ext (handles CFFI extensions)
        super().run()

        # Then compile OpenGGCM Fortran extension
        self.compile_openggcm()

    def compile_openggcm(self):
        """Compile OpenGGCM Fortran extension using f2py."""
        # Get paths
        src_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'kamodo_ccmc', 'readers', 'OpenGGCM'
        )

        # Check if Fortran source files exist
        fortran_files = ['read_b_grids.f', 'readmagfile3d.f']
        for f in fortran_files:
            if not os.path.exists(os.path.join(src_dir, f)):
                warnings.warn(f"OpenGGCM Fortran source {f} not found, skipping compilation")
                return

        # Determine output directory (where the compiled module should go)
        # During wheel build, we need to put it in build_lib
        if self.inplace:
            output_dir = src_dir
        else:
            output_dir = os.path.join(
                self.build_lib, 'kamodo_ccmc', 'readers', 'OpenGGCM'
            )
            os.makedirs(output_dir, exist_ok=True)

        # Check for existing compiled module
        existing = glob.glob(os.path.join(output_dir, 'readOpenGGCM*.so')) + \
                   glob.glob(os.path.join(output_dir, 'readOpenGGCM*.pyd'))
        if existing and not self.force:
            print(f"OpenGGCM extension already exists: {existing[0]}")
            return

        print("Compiling OpenGGCM Fortran extension with f2py...")

        # Build f2py command
        # Use numpy.f2py module directly for better compatibility
        f2py_cmd = [
            sys.executable, '-m', 'numpy.f2py',
            '-c', '-m', 'readOpenGGCM',
            'read_b_grids.f', 'readmagfile3d.f'
        ]

        try:
            # Run f2py from the source directory
            result = subprocess.run(
                f2py_cmd,
                cwd=src_dir,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                print(f"f2py stdout:\n{result.stdout}")
                print(f"f2py stderr:\n{result.stderr}")
                raise RuntimeError(f"f2py compilation failed with code {result.returncode}")

            # Find the compiled module
            compiled = glob.glob(os.path.join(src_dir, 'readOpenGGCM*.so')) + \
                       glob.glob(os.path.join(src_dir, 'readOpenGGCM*.pyd'))

            if not compiled:
                raise RuntimeError("f2py completed but no compiled module found")

            compiled_file = compiled[0]
            print(f"Compiled: {compiled_file}")

            # Move to output directory if needed
            if output_dir != src_dir:
                dest = os.path.join(output_dir, os.path.basename(compiled_file))
                shutil.move(compiled_file, dest)
                print(f"Moved to: {dest}")

            print("OpenGGCM Fortran extension compiled successfully!")

        except FileNotFoundError:
            warnings.warn(
                "f2py not found. OpenGGCM reader will not be available. "
                "Install gfortran and numpy to enable Fortran compilation."
            )
        except Exception as e:
            # Check if this is a release build (strict mode)
            if os.environ.get('KAMODO_RELEASE'):
                raise RuntimeError(f"OpenGGCM compilation failed: {e}")
            else:
                warnings.warn(
                    f"OpenGGCM Fortran compilation failed: {e}. "
                    "The OpenGGCM reader will not be available."
                )


# Register a dummy extension to ensure build_ext runs
# (CFFI handles its own extensions, but we need build_ext for Fortran)
ext_modules = [
    Extension('kamodo_ccmc.readers.OpenGGCM._openggcm_placeholder', []),
]

if __name__ == "__main__":
    setup(
        cffi_modules=cffi_modules,
        ext_modules=ext_modules,
        cmdclass={'build_ext': build_ext},
    )
