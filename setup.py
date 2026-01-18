"""Kamodo-CCMC setup script with native extension compilation.

This setup.py merges the best approaches from Claude and Codex:
- Uses CFFI's native setuptools integration (cffi_modules) for C extensions
- Compiles OpenGGCM Fortran extension via f2py subprocess (SpacePy-style)
- Bundles runtime libraries for portable wheels (SpacePy-style)
- Supports KAMODO_SKIP_NATIVE to skip native extensions
- Supports KAMODO_RELEASE for strict build error handling
"""

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import subprocess
import sys
import os
import shutil
import glob
import warnings

ROOT = os.path.abspath(os.path.dirname(__file__))
KAMODO_RELEASE = os.environ.get("KAMODO_RELEASE", "").lower() in {"1", "true", "yes"}
KAMODO_SKIP_NATIVE = os.environ.get("KAMODO_SKIP_NATIVE", "").lower() in {"1", "true", "yes"}

# Use CFFI's native setuptools integration to compile C extensions during wheel build.
# This tells setuptools to invoke each ffibuilder.compile() during the build phase.
# The format is "path/to/build_script.py:ffibuilder_object_name"
cffi_modules = [
    "kamodo_ccmc/readers/OCTREE_BLOCK_GRID/interpolate_amrdata_extension_build.py:ffibuilder",
    "kamodo_ccmc/readers/Tri2D/interpolate_tri2d_extension_build.py:ffibuilder",
]


def _get_env_with_distutils_fix():
    """Get environment with distutils compatibility for Python 3.12+.

    NumPy's f2py may need stdlib distutils on older Python versions,
    but Python 3.12+ removed it. This sets SETUPTOOLS_USE_DISTUTILS appropriately.
    """
    env = os.environ.copy()
    if "SETUPTOOLS_USE_DISTUTILS" not in env:
        if sys.version_info >= (3, 12):
            env["SETUPTOOLS_USE_DISTUTILS"] = "local"
        else:
            env["SETUPTOOLS_USE_DISTUTILS"] = "stdlib"
    return env


def _copy_macos_fortran_libs(outdir, strict=False):
    """Copy Fortran runtime libraries for portable macOS wheels (SpacePy-style).

    Args:
        outdir: Directory to copy libraries into (kamodo_ccmc package dir)
        strict: If True, raise errors; if False, warn and continue

    Returns:
        List of paths to copied library files
    """
    libnames = ["libgfortran", "libquadmath", "libgcc_s.1", "libgcc_s.1.1"]
    outlibdir = os.path.join(outdir, "libs")
    os.makedirs(outlibdir, exist_ok=True)
    outputs = []

    for lib in libnames:
        try:
            proc = subprocess.run(
                ["gfortran", f"--print-file-name={lib}.dylib"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if proc.returncode != 0:
                if strict and not lib.startswith("libgcc_s"):
                    raise RuntimeError(f"Failed locating {lib}.dylib via gfortran.")
                continue

            libpath = proc.stdout.strip()
            if not os.path.isfile(libpath):
                # libgcc_s variants may not all exist; that's OK
                if lib.startswith("libgcc_s"):
                    continue
                if strict:
                    raise RuntimeError(f"Required runtime library not found: {libpath}")
                warnings.warn(f"Skipping missing {lib}.dylib at {libpath}.")
                continue

            # Resolve symlinks to get actual file
            libpath = os.path.realpath(libpath)
            dest = os.path.join(outlibdir, os.path.basename(libpath))
            shutil.copy2(libpath, dest)
            outputs.append(dest)
            print(f"Bundled: {os.path.basename(libpath)}")

        except FileNotFoundError:
            if strict:
                raise RuntimeError("gfortran not found for library bundling.")
            warnings.warn(f"Skipping {lib}.dylib (gfortran not available).")

    return outputs


def _copy_windows_fortran_libs(outdir, strict=False):
    """Copy Fortran runtime DLLs for portable Windows wheels (SpacePy-style).

    Args:
        outdir: Directory to copy libraries into (kamodo_ccmc package dir)
        strict: If True, raise errors; if False, warn and continue

    Returns:
        List of paths to copied DLL files
    """
    outputs = []
    libneeded = ("libgfortran", "libgcc_s", "libquadmath")
    liboptional = ("libwinpthread",)
    libdir = None
    libnames = None

    # Search PATH for required DLLs
    for p in os.environ.get("PATH", "").split(os.pathsep):
        if not os.path.isdir(p):
            continue
        try:
            candidates = [
                f for f in os.listdir(p)
                if f.lower().endswith(".dll") and f.startswith(libneeded + liboptional)
            ]
        except OSError:
            continue

        # Check if all required libs are present
        if len([f for f in candidates if f.startswith(libneeded)]) == len(libneeded):
            libdir = p
            libnames = candidates
            break

    if libdir is None:
        if strict:
            raise RuntimeError("Could not locate Fortran runtime DLLs on PATH.")
        warnings.warn("Fortran runtime DLLs not found; Windows wheels may be incomplete.")
        return outputs

    outlibdir = os.path.join(outdir, "libs")
    os.makedirs(outlibdir, exist_ok=True)

    for f in libnames:
        src = os.path.join(libdir, f)
        dest = os.path.join(outlibdir, f)
        shutil.copy2(src, dest)
        outputs.append(dest)
        print(f"Bundled: {f}")

    return outputs


def _add_macos_rpath(path, strict=False):
    """Add rpath to compiled extension so it finds bundled libs.

    Args:
        path: Path to the .so file
        strict: If True, raise errors; if False, warn and continue
    """
    cmd = ["install_name_tool", "-add_rpath", "@loader_path/../libs", path]
    try:
        subprocess.check_call(cmd)
        print(f"Added rpath to: {os.path.basename(path)}")
    except (OSError, subprocess.CalledProcessError) as exc:
        if strict:
            raise RuntimeError(f"Failed adding rpath to {path}.") from exc
        warnings.warn(f"Failed adding rpath to {path}; runtime libs may not load.")


def _clean_source_artifacts():
    """Remove intermediate .o and .so files from source tree after build.

    This prevents stale build artifacts from being accidentally included
    in source distributions or contaminating future builds.
    """
    readers_dir = os.path.join(ROOT, "kamodo_ccmc", "readers")
    for root, _, files in os.walk(readers_dir):
        for filename in files:
            if filename.endswith((".o", ".so")) and not filename.startswith("_"):
                # Don't remove CFFI-generated .so files (they start with _)
                filepath = os.path.join(root, filename)
                try:
                    os.remove(filepath)
                    print(f"Cleaned: {filepath}")
                except OSError:
                    pass


class build_ext(_build_ext):
    """Custom build_ext that compiles OpenGGCM Fortran extension via f2py.

    This follows SpacePy's approach: register an Extension to trigger build_ext,
    then do the actual compilation in run() via subprocess. Also bundles runtime
    libraries for portable wheels.
    """

    def run(self):
        self._extension_outputs = []

        # Check if user wants to skip native extensions
        if KAMODO_SKIP_NATIVE:
            print("KAMODO_SKIP_NATIVE set: skipping native extension compilation.")
            return

        # First, run the standard build_ext (handles CFFI extensions)
        super().run()

        # Then compile OpenGGCM Fortran extension
        self.compile_openggcm()

        # Bundle runtime libraries for portable wheels
        self._bundle_runtime_libs()

        # Clean intermediate build artifacts from source tree
        _clean_source_artifacts()

    def compile_openggcm(self):
        """Compile OpenGGCM Fortran extension using f2py."""
        src_dir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', 'OpenGGCM')

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
            self._extension_outputs.append(existing[0])
            return

        print("Compiling OpenGGCM Fortran extension with f2py...")

        # Build f2py command
        f2py_cmd = [
            sys.executable, '-m', 'numpy.f2py',
            '-c', '-m', 'readOpenGGCM',
            'read_b_grids.f', 'readmagfile3d.f'
        ]

        try:
            # Run f2py with distutils compatibility fix
            result = subprocess.run(
                f2py_cmd,
                cwd=src_dir,
                capture_output=True,
                text=True,
                env=_get_env_with_distutils_fix()
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
                self._extension_outputs.append(dest)
            else:
                self._extension_outputs.append(compiled_file)

            print("OpenGGCM Fortran extension compiled successfully!")

        except FileNotFoundError:
            msg = ("f2py not found. OpenGGCM reader will not be available. "
                   "Install gfortran and numpy to enable Fortran compilation.")
            if KAMODO_RELEASE:
                raise RuntimeError(msg)
            warnings.warn(msg)
        except Exception as e:
            if KAMODO_RELEASE:
                raise RuntimeError(f"OpenGGCM compilation failed: {e}")
            else:
                warnings.warn(
                    f"OpenGGCM Fortran compilation failed: {e}. "
                    "The OpenGGCM reader will not be available."
                )

    def _bundle_runtime_libs(self):
        """Bundle Fortran runtime libraries for portable wheels (SpacePy-style)."""
        if self.inplace:
            outdir = os.path.join(ROOT, "kamodo_ccmc")
        else:
            outdir = os.path.join(os.path.abspath(self.build_lib), "kamodo_ccmc")

        if sys.platform == "darwin":
            libs = _copy_macos_fortran_libs(outdir, strict=KAMODO_RELEASE)
            self._extension_outputs.extend(libs)
            # Add rpath to compiled extensions so they find bundled libs
            for output in self._extension_outputs:
                if output.endswith(".so"):
                    _add_macos_rpath(output, strict=KAMODO_RELEASE)
        elif sys.platform == "win32":
            libs = _copy_windows_fortran_libs(outdir, strict=KAMODO_RELEASE)
            self._extension_outputs.extend(libs)

    def get_outputs(self):
        """Return list of output files (for proper wheel inclusion)."""
        outputs = super().get_outputs() or []
        outputs.extend(getattr(self, '_extension_outputs', []))
        return outputs


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
