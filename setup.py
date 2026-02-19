"""Kamodo-CCMC setup script with shared library compilation (SpacePy-style).

Compiles C/Fortran code directly to shared libraries (.so/.dll/.dylib) that
are loaded via ctypes at runtime. No CFFI, no f2py, no numpy build dependency.

Wheels are tagged cp310-abi3 (Stable ABI / Limited API) so a single wheel
per platform serves Python 3.10+. The authoritative abi3 declaration is
the ``options={'bdist_wheel': {'py_limited_api': 'cp310'}}`` passed to
setup() at the bottom of this file; pyproject.toml carries a secondary
belt-and-suspenders declaration.
"""

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import subprocess
import sys
import os
import shutil
import glob
import warnings

# distutils is provided by setuptools' vendored copy on Python 3.12+
# (stdlib distutils removed in 3.12; pyproject.toml ensures setuptools is present)
import distutils.ccompiler
import distutils.sysconfig

ROOT = os.path.abspath(os.path.dirname(__file__))
KAMODO_RELEASE = os.environ.get("KAMODO_RELEASE", "").lower() in {"1", "true", "yes"}
KAMODO_SKIP_NATIVE = os.environ.get("KAMODO_SKIP_NATIVE", "").lower() in {"1", "true", "yes"}


def _rewrite_macos_libgcc_dependency(libgcc_path, strict=False):
    """Rewrite hardcoded libgcc_s.1.1 path inside copied libgcc_s.1.

    Some GCC toolchains ship libgcc_s.1 with an absolute dependency on
    libgcc_s.1.1. When bundled in wheels, that absolute path can break on
    user machines. Rewrite to @loader_path-relative, matching SpacePy's guard.
    """
    basename = os.path.basename(libgcc_path)
    if sys.platform != "darwin" or not basename.startswith("libgcc_s.1"):
        return

    try:
        proc = subprocess.run(
            ["otool", "-L", libgcc_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        if strict:
            raise RuntimeError("otool not found while validating libgcc_s linkage.")
        return

    if proc.returncode != 0:
        if strict:
            raise RuntimeError("Failed to inspect libgcc_s dependencies with otool.")
        return

    hardcoded = [
        line.strip().split()[0]
        for line in proc.stdout.splitlines()
        if line.startswith("\t/") and "libgcc_s.1.1" in line
    ]

    if not hardcoded:
        return

    new_path = os.path.join("@loader_path", os.path.basename(hardcoded[0]))
    cmd = ["install_name_tool", "-change", hardcoded[0], new_path, libgcc_path]
    try:
        subprocess.check_call(cmd)
        print(f"Rewrote libgcc_s dependency: {hardcoded[0]} -> {new_path}")
    except (FileNotFoundError, subprocess.CalledProcessError):
        if strict:
            raise RuntimeError(f"Failed rewriting hardcoded libgcc path in {libgcc_path}.")


def _copy_macos_fortran_libs(outdir, strict=False):
    """Copy Fortran runtime libraries for portable macOS wheels.

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
            if lib.startswith("libgcc_s"):
                _rewrite_macos_libgcc_dependency(dest, strict=strict)
            outputs.append(dest)
            print(f"Bundled: {os.path.basename(libpath)}")

        except FileNotFoundError:
            if strict:
                raise RuntimeError("gfortran not found for library bundling.")
            warnings.warn(f"Skipping {lib}.dylib (gfortran not available).")

    return outputs


def _copy_windows_fortran_libs(outdir, strict=False):
    """Copy Fortran runtime DLLs for portable Windows wheels.

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


def _clean_source_artifacts():
    """Remove intermediate .o files from source tree after build.

    Only cleans known build intermediates from the three extension source
    directories. This prevents stale object files from contaminating
    future builds without risking deletion of legitimate library outputs.
    """
    for subdir in ['OCTREE_BLOCK_GRID', 'Tri2D', 'OpenGGCM']:
        srcdir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', subdir)
        if not os.path.isdir(srcdir):
            continue
        for filename in os.listdir(srcdir):
            if filename.endswith('.o'):
                try:
                    os.remove(os.path.join(srcdir, filename))
                except OSError:
                    pass


class build_ext(_build_ext):
    """Custom build_ext that compiles C/Fortran to shared libraries.

    C code is compiled via distutils.ccompiler (honours CC env var and
    platform defaults). Fortran is compiled via subprocess with gfortran
    (or the FC env var) since distutils has no Fortran abstraction.
    """

    def finalize_options(self):
        super().finalize_options()
        # Match SpacePy behavior: prefer MinGW toolchain by default on Windows.
        if self.compiler is None and sys.platform == "win32":
            self.compiler = "mingw32"

    @staticmethod
    def _configure_c_compiler(compiler_name):
        """Create and customize a distutils C compiler with Windows tweaks."""
        comp = distutils.ccompiler.new_compiler(compiler=compiler_name)
        if sys.platform == "win32":
            # Avoid MSVC runtime coupling when building with mingw.
            comp.dll_libraries = []
        distutils.sysconfig.customize_compiler(comp)
        return comp

    @staticmethod
    def _c_build_args(comp):
        """Return compile/link args compatible with active C compiler."""
        is_msvc = getattr(comp, "compiler_type", "") == "msvc"
        compile_args = [] if is_msvc else ["-O2"]
        link_libs = [] if is_msvc else ["m"]
        return compile_args, link_libs

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
        """Compile OCTREE_BLOCK_GRID as shared library."""
        srcdir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', 'OCTREE_BLOCK_GRID')

        if self.inplace:
            outdir = srcdir
        else:
            outdir = os.path.abspath(os.path.join(self.build_lib, 'kamodo_ccmc', 'readers', 'OCTREE_BLOCK_GRID'))
            os.makedirs(outdir, exist_ok=True)

        # Use distutils compiler abstraction (honours CC env var)
        comp = self._configure_c_compiler(self.compiler)
        libname = comp.library_filename('interpolate_amrdata', lib_type='shared')
        libpath = os.path.join(outdir, libname)

        # Check if up-to-date
        source_files = [
            'interpolate_amrdata.c', 'setup_parent.c', 'setup_octree.c',
            'interpolate_in_block.c', 'find_octree_block.c', 'find_in_block.c',
            'find_block.c', 'trace_fieldline.c', 'ctypes_wrappers.c',
        ]
        sources = [os.path.join(srcdir, f) for f in source_files]
        if os.path.exists(libpath):
            src_newest = max(os.path.getmtime(s) for s in sources)
            lib_mtime = os.path.getmtime(libpath)
            if lib_mtime > src_newest and not self.force:
                print(f"OCTREE library up-to-date: {libname}")
                return libpath

        print("Compiling OCTREE_BLOCK_GRID as shared library...")

        # Compile C sources via distutils (no hardcoded gcc)
        # -Dno_idl ensures we use the ctypes-compatible function signatures
        compile_args, link_libs = self._c_build_args(comp)
        compile_kwargs = {
            "output_dir": self.build_temp,
            "macros": [('no_idl', None)],
        }
        if compile_args:
            compile_kwargs["extra_preargs"] = compile_args
        objects = comp.compile(sources, **compile_kwargs)

        # Link to shared library; set rpath at link time on macOS
        extra_link_args = []
        if sys.platform == 'darwin':
            extra_link_args = ['-Wl,-rpath,@loader_path/../libs']

        comp.link_shared_lib(objects, 'interpolate_amrdata',
                             libraries=link_libs, output_dir=outdir,
                             extra_postargs=extra_link_args)

        print(f"Compiled: {libpath}")
        return libpath

    def compile_tri2d(self):
        """Compile Tri2D as shared library."""
        srcdir = os.path.join(ROOT, 'kamodo_ccmc', 'readers', 'Tri2D')

        if self.inplace:
            outdir = srcdir
        else:
            outdir = os.path.abspath(os.path.join(self.build_lib, 'kamodo_ccmc', 'readers', 'Tri2D'))
            os.makedirs(outdir, exist_ok=True)

        # Use distutils compiler abstraction (honours CC env var)
        comp = self._configure_c_compiler(self.compiler)
        libname = comp.library_filename('interpolate_tri2d', lib_type='shared')
        libpath = os.path.join(outdir, libname)

        # Check if up-to-date
        source_files = [
            'interpolate_tri2d_plus_1d.c', 'setup_tri.c',
            'find_tri.c', 'hunt.c', 'ctypes_wrappers.c',
        ]
        sources = [os.path.join(srcdir, f) for f in source_files]
        if os.path.exists(libpath):
            src_newest = max(os.path.getmtime(s) for s in sources)
            lib_mtime = os.path.getmtime(libpath)
            if lib_mtime > src_newest and not self.force:
                print(f"Tri2D library up-to-date: {libname}")
                return libpath

        print("Compiling Tri2D as shared library...")

        # Compile C sources via distutils (no hardcoded gcc)
        compile_args, link_libs = self._c_build_args(comp)
        compile_kwargs = {"output_dir": self.build_temp}
        if compile_args:
            compile_kwargs["extra_preargs"] = compile_args
        objects = comp.compile(sources, **compile_kwargs)

        # Link to shared library; set rpath at link time on macOS
        extra_link_args = []
        if sys.platform == 'darwin':
            extra_link_args = ['-Wl,-rpath,@loader_path/../libs']

        comp.link_shared_lib(objects, 'interpolate_tri2d',
                             libraries=link_libs, output_dir=outdir,
                             extra_postargs=extra_link_args)

        print(f"Compiled: {libpath}")
        return libpath

    def compile_openggcm(self):
        """Compile OpenGGCM Fortran as shared library.

        Fortran stays subprocess-based (distutils has no Fortran abstraction).
        """
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
            outdir = os.path.abspath(os.path.join(self.build_lib, 'kamodo_ccmc', 'readers', 'OpenGGCM'))
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

            # Link to shared library; set rpath at link time on macOS
            ldflags = ['-shared', '-fPIC']
            if sys.platform == 'darwin':
                ldflags.extend(['-mmacosx-version-min=11.0',
                                '-Wl,-rpath,@loader_path/../libs'])

            link_cmd = [fc] + ldflags + glob.glob('*.o') + ['-o', libname]
            subprocess.check_call(link_cmd)

            # Move to output directory
            if outdir != srcdir:
                shutil.move(libname, libpath)

            print(f"Compiled: {libpath}")
            return libpath

        finally:
            os.chdir(original_dir)

    def _bundle_runtime_libs(self):
        """Bundle Fortran runtime libraries for portable wheels."""
        if self.inplace:
            outdir = os.path.join(ROOT, "kamodo_ccmc")
        else:
            outdir = os.path.join(os.path.abspath(self.build_lib), "kamodo_ccmc")

        if sys.platform == "darwin":
            libs = _copy_macos_fortran_libs(outdir, strict=KAMODO_RELEASE)
            self._extension_outputs.extend(libs)
        elif sys.platform == "win32":
            libs = _copy_windows_fortran_libs(outdir, strict=KAMODO_RELEASE)
            self._extension_outputs.extend(libs)

    def get_outputs(self):
        """Return list of output files (for proper wheel inclusion)."""
        outputs = list(getattr(self, '_extension_outputs', []))
        outputs.extend(super().get_outputs())
        return outputs


# Dummy extensions to trigger build_ext
# The actual compilation is handled by our custom build_ext.run() method.
# Without these, setuptools skips build_ext entirely.
ext_modules = [
    Extension('kamodo_ccmc.readers.OCTREE_BLOCK_GRID.interpolate_amrdata', []),
    Extension('kamodo_ccmc.readers.Tri2D.interpolate_tri2d', []),
    Extension('kamodo_ccmc.readers.OpenGGCM.readOpenGGCM', []),
]


if __name__ == "__main__":
    setup(
        cmdclass={'build_ext': build_ext},
        ext_modules=ext_modules,
        options={'bdist_wheel': {'py_limited_api': 'cp310'}},
    )
