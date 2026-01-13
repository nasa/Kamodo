from setuptools import setup

# Import custom build commands
try:
    from kamodo_ccmc.build_tools import build_ext
    cmdclass = {'build_ext': build_ext}
except ImportError:
    # Fallback if build_tools can't be imported (shouldn't happen in normal builds)
    import warnings
    warnings.warn("Could not import kamodo_ccmc.build_tools, extensions won't be auto-compiled")
    cmdclass = {}

if __name__ == "__main__":
    setup(cmdclass=cmdclass)
