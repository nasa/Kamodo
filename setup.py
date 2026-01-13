from setuptools import setup

# Use CFFI's native setuptools integration to compile C extensions during wheel build.
# This tells setuptools to invoke each ffibuilder.compile() during the build phase.
# The format is "path/to/build_script.py:ffibuilder_object_name"
cffi_modules = [
    "kamodo_ccmc/readers/OCTREE_BLOCK_GRID/interpolate_amrdata_extension_build.py:ffibuilder",
    "kamodo_ccmc/readers/Tri2D/interpolate_tri2d_extension_build.py:ffibuilder",
]

if __name__ == "__main__":
    setup(cffi_modules=cffi_modules)
