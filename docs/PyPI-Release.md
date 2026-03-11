# PyPI Release Checklist

This checklist is modeled after SpacePy's `DISTRIBUTE` workflow: build a clean sdist on
Unix, then build platform wheels, test, and upload. Kamodo includes native
extensions (ctypes-loaded shared libraries compiled from C/Fortran sources), so wheels are the primary delivery
mechanism for "pip install just works."

## Environment Variables

- **`KAMODO_RELEASE=1`**: Strict mode - build fails if any native extension fails to compile.
  Always set this when building release wheels.
- **`KAMODO_SKIP_NATIVE=1`**: Skip native extension compilation entirely (for debugging only).

## 1) Prepare the Release

- [ ] Update version in `setup.cfg` and `kamodo_ccmc/__init__.py`
- [ ] Update CHANGELOG or release notes
- [ ] Ensure all tests pass on target Python versions (3.10, 3.11, 3.12, 3.13)
- [ ] Tag the release commit: `git tag v25.X.Y`

## 2) Build the Source Distribution (Unix)

Use a clean build environment:

```bash
# Clean any previous builds
rm -rf build/ dist/ *.egg-info kamodo_ccmc.egg-info

# Build sdist
python -m build -s
```

Verify the sdist contains native sources:
- C files from `kamodo_ccmc/readers/OCTREE_BLOCK_GRID/` and `kamodo_ccmc/readers/Tri2D/`
- Fortran files from `kamodo_ccmc/readers/OpenGGCM/`

## 3) Build Wheels (Multiple Platforms)

SpacePy builds wheels on Linux, macOS, and Windows so users don't need local
compilers. Kamodo follows the same approach using `cibuildwheel`.

### Option A: GitHub Actions (Recommended)

Push a version tag to trigger the wheel build workflow:

```bash
git tag v25.X.Y
git push origin v25.X.Y
```

Or manually trigger via GitHub Actions UI (workflow_dispatch).

### Option B: Local Build with cibuildwheel

```bash
pip install cibuildwheel

# Linux (in Docker)
KAMODO_RELEASE=1 python -m cibuildwheel --platform linux --output-dir wheelhouse

# macOS
KAMODO_RELEASE=1 python -m cibuildwheel --platform macos --output-dir wheelhouse

# Windows (requires MSYS2 gfortran in PATH)
set KAMODO_RELEASE=1
python -m cibuildwheel --platform windows --output-dir wheelhouse
```

### Platform Requirements

| Platform | Requirements |
|----------|-------------|
| Linux | gcc, gfortran (installed in manylinux container) |
| macOS | Xcode Command Line Tools, `brew install gcc` |
| Windows | MSYS2 with mingw-w64-x86_64-gcc-fortran |

## 4) Test the Wheels

For each platform, create a fresh environment and test:

```bash
# Create fresh environment
python -m venv test_env
source test_env/bin/activate  # or test_env\Scripts\activate on Windows

# Install wheel
pip install wheelhouse/kamodo_ccmc-*.whl

# Test imports
python -c "
import kamodo_ccmc
print(f'Version: {kamodo_ccmc.__version__}')

# Test native shared library availability (ctypes loaders)
from kamodo_ccmc.readers.OCTREE_BLOCK_GRID import OCTREE_AVAILABLE
from kamodo_ccmc.readers.Tri2D import TRI2D_AVAILABLE
from kamodo_ccmc.readers.OpenGGCM import OPENGGCM_AVAILABLE

print(f'OCTREE available: {OCTREE_AVAILABLE}')
print(f'Tri2D available: {TRI2D_AVAILABLE}')
print(f'OpenGGCM available: {OPENGGCM_AVAILABLE}')

if not all([OCTREE_AVAILABLE, TRI2D_AVAILABLE, OPENGGCM_AVAILABLE]):
    raise SystemExit('ERROR: One or more native readers are unavailable')

# Test reader modules import against installed wheel
from kamodo_ccmc.readers import swmfgm_4D, gameragm_4D
print('Reader modules: OK')

print('All extensions loaded successfully!')
"
```

## 5) Upload to PyPI

### Test PyPI First

```bash
pip install twine

# Upload to TestPyPI
twine upload --repository testpypi dist/* wheelhouse/*

# Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ kamodo-ccmc
```

### Production PyPI

After verifying TestPyPI installation works:

```bash
# Upload to production PyPI
twine upload dist/* wheelhouse/*
```

Or use GitHub Actions with trusted publishing (OIDC) - no API token needed.

## 6) Post-Release

- [ ] Create GitHub Release with changelog
- [ ] Verify `pip install kamodo-ccmc` works in fresh environments
- [ ] Announce on PyHC mailing list / Slack if appropriate

## Troubleshooting

### Extension compilation fails

- Ensure gcc/gfortran are installed and in PATH
- For source installs, use `pip install . -v` and inspect the native build logs
- Set `KAMODO_SKIP_NATIVE=1` to skip extensions (for debugging)

### Wheel missing runtime libraries

- On macOS: Ensure `brew install gcc` was run before building
- On Windows: Ensure MSYS2 gfortran DLLs are in PATH
- The build process should bundle libgfortran, libquadmath, etc.

### Import fails on installed wheel

- Check that the `libs/` directory exists in the installed package
- On Windows, verify DLLs were bundled and PATH is set correctly

## References

- [SpacePy DISTRIBUTE guide](https://github.com/spacepy/spacepy/blob/main/DISTRIBUTE)
- [cibuildwheel documentation](https://cibuildwheel.readthedocs.io/)
- [PyPA publishing guide](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/)
