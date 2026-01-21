"""Kamodo-CCMC setup script with shared library compilation.

This setup.py follows the SpacePy approach (PR #749) to compile native extensions:
- Compiles C/Fortran code directly to shared libraries (.so/.dll/.dylib)
- No CFFI, no f2py, no numpy build dependency
- Python-version-independent wheels (cp36-abi3)
- Libraries loaded via ctypes at runtime
- Portable wheels with bundled runtime libraries
"""

from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext
import subprocess
import sys
import os
import shutil
import glob
import warnings
import distutils.ccompiler
import distutils.sysconfig

ROOT = os.path.abspath(os.path.dirname(__file__))
KAMODO_RELEASE = os.environ.get("KAMODO_RELEASE", "").lower() in {"1", "true", "yes"}
KAMODO_SKIP_NATIVE = os.environ.get("KAMODO_SKIP_NATIVE", "").lower() in {"1", "true", "yes"}


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
        path: Path to the .so/.dylib file
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
            if filename.endswith((".o", ".so", ".dylib", ".dll")) and not filename.startswith("lib"):
                filepath = os.path.join(root, filename)
                try:
                    os.remove(filepath)
                    print(f"Cleaned: {filepath}")
                except OSError:
                    pass


class build_ext(_build_ext):
    """Custom build_ext that compiles extensions to shared libraries (SpacePy-style).

    Unlike CFFI/f2py approach, this compiles C/Fortran directly with gcc/gfortran
    to shared libraries (.so/.dll/.dylib) that are loaded via ctypes at runtime.
    """

    def run(self):
        self._extension_outputs = []

        # Check if user wants to skip native extensions
        if KAMODO_SKIP_NATIVE:
            print("KAMODO_SKIP_NATIVE set: skipping native extension compilation.")
            return

        # Compile OCTREE_BLOCK_GRID C extension
        try:
            octree_lib = self.compile_octree()
            if octree_lib:
                self._extension_outputs.append(octree_lib)
        except Exception as e:
            if KAMODO_RELEASE:
                raise
            warnings.warn(f"OCTREE compilation failed: {e}. SWMF-GM reader will not be available.")

        # Compile Tri2D C extension
        try:
            tri2d_lib = self.compile_tri2d()
            if tri2d_lib:
                self._extension_outputs.append(tri2d_lib)
        except Exception as e:
            if KAMODO_RELEASE:
                raise
            warnings.warn(f"Tri2D compilation failed: {e}. GAMERA-GM reader will not be available.")

        # Compile OpenGGCM Fortran extension
        try:
            openggcm_lib = self.compile_openggcm()
            if openggcm_lib:
                self._extension_outputs.append(openggcm_lib)
        except Exception as e:
            if KAMODO_RELEASE:
                raise
            warnings.warn(f"OpenGGCM compilation failed: {e}. OpenGGCM reader will not be available.")

        # Bundle runtime libraries for portable wheels
        self._bundle_runtime_libs()

        # Clean intermediate build artifacts from source tree
        _clean_source_artifacts()

    def compile_octree(self):
        """Compile OCTREE_BLOCK_GRID as shared library (SpacePy-style)."""
        srcdir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', 'OCTREE_BLOCK_GRID')

        if self.inplace:
            outdir = srcdir
        else:
            outdir = os.path.join(self.build_lib, 'kamodo_ccmc', 'readers', 'OCTREE_BLOCK_GRID')
            os.makedirs(outdir, exist_ok=True)

        # Determine shared library filename
        ccomp = distutils.ccompiler.new_compiler(compiler=self.compiler)
        distutils.sysconfig.customize_compiler(ccomp)
        libname = ccomp.library_filename('interpolate_amrdata', lib_type='shared')
        libpath = os.path.join(outdir, libname)

        # Check if up-to-date
        sources = glob.glob(os.path.join(srcdir, '*.c'))
        if os.path.exists(libpath):
            src_newest = max(os.path.getmtime(s) for s in sources)
            lib_mtime = os.path.getmtime(libpath)
            if lib_mtime > src_newest and not self.force:
                print(f"OCTREE library up-to-date: {libname}")
                return libpath

        print("Compiling OCTREE_BLOCK_GRID as shared library...")

        # Compile .c files to .o files
        original_dir = os.getcwd()
        try:
            os.chdir(srcdir)

            c_files = [
                'interpolate_amrdata.c',
                'setup_parent.c',
                'setup_octree.c',
                'interpolate_in_block.c',
                'find_octree_block.c',
                'find_in_block.c',
                'find_block.c',
                'trace_fieldline.c',
            ]

            gcc_cmd = ['gcc', '-c', '-O2', '-fPIC'] + c_files
            subprocess.check_call(gcc_cmd)

            # Link to shared library
            link_cmd = ['gcc', '-shared', '-fPIC'] + glob.glob('*.o') + ['-o', libname, '-lm']
            subprocess.check_call(link_cmd)

            # Move to output directory
            if outdir != srcdir:
                shutil.move(libname, libpath)

            # Add rpath on macOS
            if sys.platform == 'darwin':
                _add_macos_rpath(libpath, strict=KAMODO_RELEASE)

            print(f"Compiled: {libpath}")
            return libpath

        finally:
            os.chdir(original_dir)

    def compile_tri2d(self):
        """Compile Tri2D as shared library (SpacePy-style)."""
        srcdir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', 'Tri2D')

        if self.inplace:
            outdir = srcdir
        else:
            outdir = os.path.join(self.build_lib, 'kamodo_ccmc', 'readers', 'Tri2D')
            os.makedirs(outdir, exist_ok=True)

        # Determine shared library filename
        ccomp = distutils.ccompiler.new_compiler(compiler=self.compiler)
        distutils.sysconfig.customize_compiler(ccomp)
        libname = ccomp.library_filename('interpolate_tri2d', lib_type='shared')
        libpath = os.path.join(outdir, libname)

        # Check if up-to-date
        sources = glob.glob(os.path.join(srcdir, '*.c'))
        if os.path.exists(libpath):
            src_newest = max(os.path.getmtime(s) for s in sources)
            lib_mtime = os.path.getmtime(libpath)
            if lib_mtime > src_newest and not self.force:
                print(f"Tri2D library up-to-date: {libname}")
                return libpath

        print("Compiling Tri2D as shared library...")

        # Compile .c files to .o files
        original_dir = os.getcwd()
        try:
            os.chdir(srcdir)

            c_files = [
                'interpolate_tri2d_plus_1d.c',
                'setup_tri.c',
                'find_tri.c',
                'hunt.c',
            ]

            gcc_cmd = ['gcc', '-c', '-O2', '-fPIC'] + c_files
            subprocess.check_call(gcc_cmd)

            # Link to shared library
            link_cmd = ['gcc', '-shared', '-fPIC'] + glob.glob('*.o') + ['-o', libname, '-lm']
            subprocess.check_call(link_cmd)

            # Move to output directory
            if outdir != srcdir:
                shutil.move(libname, libpath)

            # Add rpath on macOS
            if sys.platform == 'darwin':
                _add_macos_rpath(libpath, strict=KAMODO_RELEASE)

            print(f"Compiled: {libpath}")
            return libpath

        finally:
            os.chdir(original_dir)

    def compile_openggcm(self):
        """Compile OpenGGCM Fortran as shared library (SpacePy-style)."""
        srcdir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', 'OpenGGCM')

        # Check if Fortran source files exist
        fortran_files = ['read_b_grids.f', 'readmagfile3d.f']
        for f in fortran_files:
            if not os.path.exists(os.path.join(srcdir, f)):
                warnings.warn(f"OpenGGCM Fortran source {f} not found, skipping compilation")
                return None

        if self.inplace:
            outdir = srcdir
        else:
            outdir = os.path.join(self.build_lib, 'kamodo_ccmc', 'readers', 'OpenGGCM')
            os.makedirs(outdir, exist_ok=True)

        # Determine shared library filename
        libname = {'darwin': 'libreadOpenGGCM.dylib',
                   'win32': 'readOpenGGCM.dll',
                   'linux': 'libreadOpenGGCM.so'}.get(sys.platform, 'libreadOpenGGCM.so')
        libpath = os.path.join(outdir, libname)

        # Check if up-to-date
        sources = [os.path.join(srcdir, f) for f in fortran_files]
        if os.path.exists(libpath):
            src_newest = max(os.path.getmtime(s) for s in sources)
            lib_mtime = os.path.getmtime(libpath)
            if lib_mtime > src_newest and not self.force:
                print(f"OpenGGCM library up-to-date: {libname}")
                return libpath

        print("Compiling OpenGGCM Fortran as shared library...")

        # Fortran compiler from environment or default
        fc = os.environ.get('FC', 'gfortran.exe' if sys.platform == 'win32' else 'gfortran')

        # Compile Fortran to .o files
        original_dir = os.getcwd()
        try:
            os.chdir(srcdir)

            compile_cmd = [fc, '-c', '-w', '-O2', '-fPIC', '-ffixed-line-length-none', '-std=legacy'] + fortran_files
            subprocess.check_call(compile_cmd)

            # Link to shared library
            ldflags = ['-shared', '-fPIC']
            if sys.platform == 'darwin':
                ldflags.extend(['-mmacosx-version-min=11.0'])

            link_cmd = [fc] + ldflags + glob.glob('*.o') + ['-o', libname]
            subprocess.check_call(link_cmd)

            # Move to output directory
            if outdir != srcdir:
                shutil.move(libname, libpath)

            # Add rpath on macOS
            if sys.platform == 'darwin':
                _add_macos_rpath(libpath, strict=KAMODO_RELEASE)

            print(f"Compiled: {libpath}")
            return libpath

        finally:
            os.chdir(original_dir)

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
                if output.endswith((".so", ".dylib")):
                    # Only add rpath if not already added
                    if not any(output.endswith(lib.split('/')[-1]) for lib in libs):
                        _add_macos_rpath(output, strict=KAMODO_RELEASE)
        elif sys.platform == "win32":
            libs = _copy_windows_fortran_libs(outdir, strict=KAMODO_RELEASE)
            self._extension_outputs.extend(libs)

    def get_outputs(self):
        """Return list of output files (for proper wheel inclusion)."""
        outputs = getattr(self, '_extension_outputs', [])
        return outputs


if __name__ == "__main__":
    setup(
        cmdclass={'build_ext': build_ext},
    )
