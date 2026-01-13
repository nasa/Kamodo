from setuptools import setup
import sys
import os

# Add the package directory to the path so we can import build_tools
# This is needed because during pip install, the package isn't installed yet
pkg_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, pkg_dir)

# Import custom build commands
try:
    from kamodo_ccmc.build_tools import build_ext, install, develop
    cmdclass = {
        'build_ext': build_ext,
        'install': install,
        'develop': develop,
    }
except ImportError as e:
    # Fallback if build_tools can't be imported
    import warnings
    warnings.warn(f"Could not import kamodo_ccmc.build_tools ({e}), extensions won't be auto-compiled")
    cmdclass = {}
finally:
    # Clean up sys.path
    if pkg_dir in sys.path:
        sys.path.remove(pkg_dir)

if __name__ == "__main__":
    setup(cmdclass=cmdclass)
