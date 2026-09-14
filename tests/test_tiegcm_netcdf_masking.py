"""Small local netCDF fixtures for TIEGCM's standard and density readers."""
from types import SimpleNamespace
from unittest.mock import patch

from netCDF4 import Dataset
import numpy as np
import pytest

from kamodo_ccmc.readers import reader_utilities as RU
from kamodo_ccmc.readers.tiegcm_4D import MODEL, model_varnames


@pytest.mark.parametrize('file_format', ['NETCDF3_CLASSIC', 'NETCDF4'])
@pytest.mark.parametrize('gvar', ['TEC', 'TN', 'DEN'])
@pytest.mark.parametrize('file_count,time_count', [(1, 2), (2, 2), (1, 1)])
def test_large_missing_values(tmp_path, file_format, gvar, file_count, time_count):
    """Use both real backends, including the next-file slice and log interpolation."""
    fill = np.float32(1e36)
    files = []
    for index in range(file_count):
        filename = str(tmp_path / f's{index + 1:03d}.nc')
        files.append(filename)
        with Dataset(filename, 'w', format=file_format) as dataset:
            for name, size in [('time', time_count), ('level', 2), ('lat', 3), ('lon', 4)]:
                dataset.createDimension(name, size)
            dimensions = ('time', 'lat', 'lon') if gvar == 'TEC' else ('time', 'level', 'lat', 'lon')
            variable = dataset.createVariable(gvar, 'f4', dimensions, fill_value=fill)
            variable.missing_value = fill
            shape = tuple(len(dataset.dimensions[name]) for name in dimensions)
            data = np.full(shape, 10.0 * (index + 1), dtype=np.float32)
            # Mask the last longitude, leaving the first longitude usable for interpolation.
            data[..., -1] = fill
            variable[:] = data

    reader_type = MODEL()
    varname, _, _, _, _, _, units = model_varnames[gvar]
    reader = SimpleNamespace(
        variables={varname: {'units': units, 'data': gvar}},
        pattern_files={'s': files},
        times={'s': {
            'all': np.arange(file_count * time_count, dtype=float),
            'start': np.arange(file_count) * time_count,
            'end': np.arange(file_count) * time_count + time_count - 1,
            'start_index': list(range(0, file_count * time_count + 1, time_count)),
        }},
        _lon=np.array([-180., -90., 0., 90., 180.]),
        _lat=np.array([-90., -45., 0., 45., 90.]),
        _ilev=np.array([1., 2.]),
        _ilev1=np.array([1., 2.]),
        ilev_sub=False,
        total_ilev=[],
    )
    reader.wrap_3Dlatlon = lambda name, data: reader_type.wrap_3Dlatlon(reader, name, data)
    reader.wrap_4Dlatlon = lambda name, data: reader_type.wrap_4Dlatlon(reader, name, data)

    # Capture the loader from the real registration method without requiring global
    # model metadata or pressure-to-height preprocessing. File IO and wrapping are real.
    with patch.object(RU, 'Functionalize_Dataset', return_value=reader) as register:
        reader_type.register_variable(reader, varname, False)
    load = register.call_args.kwargs['func']
    result = load(0)

    if gvar == 'DEN':
        # The density loader creates a real SciPy log-space interpolator.
        point = [-180., 0., 1.5]
        if time_count > 1:
            point.insert(0, 0.5)
        assert result([point])[0] == pytest.approx(10.0)
        if file_count > 1:
            assert result([[2., -180., 0., 1.5]])[0] == pytest.approx(20.0)
        missing = [90., 0., 1.5]
        if time_count > 1:
            missing.insert(0, 0.5)
        assert np.isnan(result([missing])[0])
    else:
        assert np.isnan(result).any()
        assert np.nanmax(result) == 10.0 * file_count
        assert not np.any(result >= fill)
        if time_count > 1:
            assert result.shape[0] == time_count + (file_count > 1)
            if file_count > 1:
                # This is the first sample appended from the next file.
                assert np.nanmax(result[-1]) == 20.0
                assert np.isnan(result[-1]).any()
