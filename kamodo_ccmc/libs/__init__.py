# Package for bundled Fortran runtime libraries.
#
# This directory contains runtime libraries (libgfortran, libquadmath, etc.)
# bundled during wheel builds for portable installation without requiring
# users to have gfortran installed.
#
# On macOS: .dylib files with rpath set on compiled extensions
# On Windows: .dll files added to PATH and os.add_dll_directory()
