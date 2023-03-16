# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 10:18:13 2021

@author: rringuet
version 3: dimension agnostic
version 3a: separated offsets, adapted to new SF structure
version 3b: update to new version of SF code, enable time reconstruction
    options.
version 3c: fine-tune code based on testing, add nan feature
version 3d: add averaging of model slices over two dimensions not reconstructed
version 4: implement coordinate-dependent plot labeling
Public version 0: fine-tuned orbit slicing to produce correct answer by
    subtracting longitude per time (15 deg/hr) for imaginary satellites in
    model orbit slicing to avoid shifting the local time of the new
    'satellites'.
"""

import numpy as np
import time as ti
from csv import DictReader
from datetime import datetime, timezone
from kamodo import Kamodo, kamodofy, gridify
import kamodo_ccmc.flythrough.SatelliteFlythrough as SF
from scipy.interpolate import RegularGridInterpolator
import kamodo_ccmc.flythrough.model_wrapper as MW


def read_GDC_sattraj(filename):
    # read data from Bob's file into dictionary
    file_dict = {}
    sat_file = DictReader(open(filename, 'r'), delimiter='\t')
    for row in sat_file:
        for column, value in row.items():
            file_dict.setdefault(column, []).append(float(value))
    for key in file_dict.keys():
        if key in ['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second',
                   'Millisecond']:
            file_dict[key] = np.array(file_dict[key], dtype=int)
        else:
            file_dict[key] = np.array(file_dict[key], dtype=float)

    # convert time data into utc timestamps (millisec all equal 0.0, ignore)
    file_dict['UTCtimestamps'] = time_to_utcts(
        file_dict['Year'], file_dict['Month'], file_dict['Day'],
        file_dict['Hour'], file_dict['Minute'], file_dict['Second'])
    return file_dict


@np.vectorize
def time_to_utcts(year, month, day, hour, minute, second):
    '''converts given values into a utc timestamp'''
    dt = datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc)
    return datetime.timestamp(dt)


def ts_to_dt(time_val):
    '''Convert UTC timestamp to a UTC datetime object.'''
    return datetime.utcfromtimestamp(time_val).replace(tzinfo=timezone.utc)


@np.vectorize
def ts_to_hours(time_val):
    '''return the utc hour from a utc timestamp.'''
    return datetime.utcfromtimestamp(time_val).hour


def ts_to_LT(time_val, longitude):
    '''Converts the utc timestamp into a local time using longitude.'''
    # convert to local time
    hours = ts_to_hours(time_val)
    LT = np.array(np.floor(longitude/15.), dtype=int) + hours

    # correct local times to be from 0 to 23
    neg_idx = np.where(LT < 0)[0]
    LT[neg_idx] += 24
    high_idx = np.where(LT > 23)[0]
    LT[high_idx] -= 24
    return LT


class Constellation(Kamodo):
    '''Uses shifted copies of the satellite trajectory to reconstruct
    the model data for the variable name given. See
    https://doi.org/10.22541/essoar.167214257.73153757/v1 for more info.

    - model: 'CTIPe', 'IRI', 'GITM', 'SWMF_IE', 'TIEGCM', etc
        (See MW.Choose_Model() output.)
    - variable_name = choose from list of standardized variable names from
        model chosen
    - file_dir = file path to model data
    - sat_time: a numpy array of the utc timestamps
    - c1, c2, c3: numpy arrays of the positions correlating to the utc
        timestamps (c1, c2, c3) should be (x,y,z) in R_E for cartesian
        coordinates, and (lon, lat, radius (R_E) or altitude (km)) for
        spherical coordinates.
    - coord_type: one of 'GDZ', 'GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG', 'SPH',
        'RLL' or those supported by AstroPy (See flythrough.utils.CoordConv)
    - coord_grid: either 'car' or 'sph'. Note that not all combinations
        make sense (e.g. 'SPH' and 'car') and are not allowed.
        All cartesian coordinate systems expect x, y, z values
        in R_E. All times should be in UTC timestamps.
        Note: If the input coordinates are spherical and you include offsets in
            latitude, make sure the offsets in latitude do not result in
            latitudes outside of the typical (-90,90) range in degrees. The
            program cannot compensate and will error.
    - recon_option = string representing reconstruction and comparison method
        Choose from 'UnMod_AvgSlice','UnMod_AvgDSlice','AvgMod_AvgSlice',
        'AvgMod_AvgDSlice', 'UnMod_OrbitSliceD', 'UnMod_OrbitSliceN',
        'AvgMod_OrbitSliceD', and 'AvgMod_OrbitSliceN'
        - UnMod means the satellite trajectory is not modified
        - AvgMod means the dimensions not reconstructed are replaced with
          either the average value of the original and offset satellite
          trajectories or the value given in the xx_avg keywords.
        - AvgSlice means recon is compared with a slice at the average
          dimensions not reconstructed, again using either the average or the
          given value in the xx_avg keywords.
        - AvgDSlice means recon is compared with an average of all the slices
          in the range of the dimensions not reconstructed. Weighted averages
          are not supported.
        - OrbitSlice means the model data is retrieved along all trajectories
          in the given constellation instead of along a given dimensional
          average.
          The OrbitSlice option is what an inifite number of satellites would
          see if positioned along the trajectory with time offsets t+dt, where
          dt is the difference in time between satellites along the orbit. For
          computation purposes, the user may enter the value of dt in seconds
          in the input line or use the default of 60 seconds. The trajectories
          are split into orbits, defined as where the latitude differential
          changes from positive to negative, and then extracted from the model
          for the calculated time grid. The 'D' and 'N' variations are for
          retrieving only the day or night portions of the orbit rather than
          the entire orbit to avoid averaging data of different natures.
          Whether a given time+spatial coordinate is on the day or night side
          is determined by converting the time+spatial coordinates of the orbit
          slices into GSE spherical coordinates and filtering based on
          longitude. The longitude logic compensates to keep equatorial
          crossing in GSE at the same MLT as the original set.
          NOTE: Orbital slicing may not produce the desired result in cartesian
          coordinates.
    - recon_dims = string representing choice of dimensions to use in
        reconstruction. 'tc1','tc2','tc3','c1c2','c1c3','c2c3' are the options.
        The first name is the x-axis, and the second name is the y-axis.
    - time_offsets: a list of time offsets in seconds to be applied to the
        trajectories (in s)
    - c1_offsets: a list of offsets to be applied to the c1 coordinate in the
        trajectory in the same units indicated by the coordinate system options
    - c2_offsets: a list of offsets to be applied to the c2 coordinate in the
        trajectory in the same units indicated by the coordinate system options
    - c3_offsets: a list of offsets to be applied to the c3 coordinate in the
        trajectory in the same units indicated by the coordinate system options
    - the offset lists entered MUST BE THE SAME LENGTH or be equal to [0.], and
        must be in the same coordinate system and units as the input
        coordinates.
    - t_avg: the average time in UTC seconds used for when time is not a
        reconstructed dimension and the recon_option chosen requires averaging.
        Default is the average time of the entire set of input and offset
        trajectories.
    - c1_avg: the average c1 value used for when c1 is not a reconstructed
        dimension and the recon_option chosen requires averaging. Default is
        the average c1 value of the entire set of input and offset
        trajectories.
    - c2_avg: the average c2 value used for when c2 is not a reconstructed
        dimension and the recon_option chosen requires averaging. Default is
        the average c2 value of the entire set of input and offset
        trajectories.
    - c3_avg: the average c3 value used for when c3 is not a reconstructed
        dimension and the recon_option chosen requires averaging. Default is
        the average c3 value of the entire set of input and offset
        trajectories.
    - dx = grid resolution for the x axis in the same units as x, default of 2.
    - dy = grid resolution for the y axis in the same units as y, default of 2.
    - d1 = grid resolution for the first non-reconstructed dimension. Only used
        if *_AvgDSlice options are chosen and if slicing in local time.
        Default=0.
    - d2 = grid resolution for the second non-reconstructed dimension. Only
        used if *_AvgDSlice options are chosen. Default=0.
    - dt = time resolution of sampling for orbit slicing option. Default=60 s.
    - run_option = one of 'all' or 'flythrough'. Default is all.
    Note: The dx and dy values are assumed to be in the units of the associated
        coordinate indicated by the coordinate options (s, R_E, deg, or km).
    Note: When reconstructing in two spatial cartesian coordinates
        (e.g. xy, yz, or xz), filter the input satellite trajectory to only
        include times where the third spatial dimension inputs are either
        positive OR negative values to avoid a combined reconstruction of both
        hemispheres in the same plots. For example, if you wish to perform a
        reconstruction of the x and y dimensions in spatial coordinates for
        the North Pole, then filter the satellite trajectory to only include
        times where the z dimension is positive. Similarly, reconstructions of
        the South Pole can be performed by simply filtering out the times where
        the z dimension is positive. This filtering is not needed when time is
        one of the reconstructed dimensions, but it can be done if desired.
        Example:
            import numpy as np
            pos_idx = np.where(results['c3']>0.)[0]  # North Pole is z>0.
            time = results['utc_time'][pos_idx]
            c1 = results['c1'][pos_idx]
            c2 = results['c2'][pos_idx]
            c3 = results['c3'][pos_idx]
    Cartesian note: For cartesian reconstructions of nearly spherical
        trajectories, 'Unmod_AvgDSlice' option is recommended to avoid
        reconstructing a ring of data around the Earth. Reconstructions will be
        better achieved in spherical coordinates for these scenarios. To easily
        switch between coordinate systems, from GSE cartesian to GEO spherical,
        for example:
            from kamodo_ccmc.flythrough.utils import ConvertCoord
            lon, lat, radius, units = ConvertCoord(utc_time, x, y, z,
                          'GSE', 'car', 'GEO', 'sph')
        See the relevant flythrough example notebook for more details on this
        function.
    '''

    def __init__(self, model, variable_name, file_dir, sat_time, c1,
                 c2, c3, coord_type, coord_grid, recon_option, recon_dims,
                 time_offsets=[0.], c1_offsets=[0.], c2_offsets=[0.],
                 c3_offsets=[0.], dx=2., dy=2., d1=0., d2=0., LT='',
                 t_avg='', c1_avg='', c2_avg='', c3_avg='', dt=60.,
                 run_option='all', **kwargs):
        super(Constellation, self).__init__(**kwargs)
        t_start = ti.perf_counter()  # Begin execution timer

        # store input values
        self.model = model
        self.variable_name = variable_name
        self.file_dir = file_dir  # file path to model data
        self.dx, self.dy = dx, dy  # grid resolution for x and y axes
        self.d1, self.d2 = d1, d2  # resolution for non-reconstructed dims
        # save input satellite trajectory time (utc timestamps)
        self.sat_time = sat_time
        self.c1, self.c2, self.c3 = c1, c2, c3  # save input spatial trajectory
        self.coord_type = coord_type
        self.coord_grid = coord_grid  # 'car' or 'sph'
        self.recon_option = recon_option
        self.recon_dims = recon_dims  # which two dimensions to reconstruct
        # datetime object for first utc timestamp
        self.sat_datetime = ts_to_dt(sat_time[0])
        self.t_avg, self.c1_avg, self.c2_avg, self.c3_avg = t_avg, c1_avg,\
            c2_avg, c3_avg
        self.dt = dt  # used for orbit slicing
        if 't' in recon_dims:
            self.dx = dx/3600.  # change dx from sec to hours
        self.run_option = run_option  # either only flythrough or all

        # check inputs for validity
        if recon_option not in ['UnMod_AvgSlice', 'UnMod_AvgDSlice',
                                'AvgMod_AvgSlice', 'AvgMod_AvgDSlice',
                                'UnMod_OrbitSliceD', 'UnMod_OrbitSliceN',
                                'AvgMod_OrbitSliceD', 'AvgMod_OrbitSlice_N']:
            print("Possible values for recon_option are 'UnMod_AvgSlice', " +
                  "'UnMod_AvgDSlice', 'AvgMod_AvgSlice', 'AvgMod_AvgDSlice'," +
                  " 'UnMod_OrbitSliceD', 'UnMod_OrbitSliceN', " +
                  "'AvgMod_OrbitSliceD', and 'AvgMod_OrbitSlice_N'.")
            raise AttributeError('Value for recon_option not recognized.')
        if recon_dims not in ['tc1', 'tc2', 'tc3', 'c1c2', 'c1c3', 'c2c3']:
            print("Possible values for recon_dims are 'tc1', 'tc2', 'tc3', " +
                  "'c1c2', 'c1c3', or 'c2c3'.")
            raise AttributeError('Value for recon_dims not recognized.')
        if 'AvgDSlice' in recon_option and (d1 == 0. or d2 == 0.):
            raise AttributeError('You must specify d1 and d2 for ' +
                                 f'{self.recon_option}.')

        # sometimes see [0.,0.,0.] as an offset for no input values (?!).
        # seems to be a memory random error. Doesn't happen all the time.
        if sum(time_offsets) == 0:
            time_offsets = [0.]
        if sum(c1_offsets) == 0:
            c1_offsets = [0.]
        if sum(c2_offsets) == 0:
            c2_offsets = [0.]
        if sum(c3_offsets) == 0:
            c3_offsets = [0.]

        # check that offset lists are the same length. append zeros if not
        # print(time_offsets, c1_offsets, c2_offsets, c3_offsets)
        max_length = max([len(item) for item in [time_offsets, c1_offsets,
                                                 c2_offsets, c3_offsets]])
        for item in [time_offsets, c1_offsets, c2_offsets, c3_offsets]:
            while len(item) < max_length:
                item.append(0.)
        self.time_offsets = time_offsets
        self.c1_offsets = c1_offsets
        self.c2_offsets = c2_offsets
        self.c3_offsets = c3_offsets

        # get variable units from satellite flythrough software
        self.variable_units = MW.Var_units(model,
                                           [variable_name])[variable_name]
        self.variables = {}   # initialize kamodo's dictionary

        # BEGIN ANALYSIS
        if self.run_option == 'all':
            # fly through model data, accounting for shifts and wrapping of lon
            results = self.flythrough()  # full dictionary from flythrough
            # sort data onto a grid, assign interpolator
            self.recon_grid = self.grid_data(results)
            self.create_2Dinterpolator(self.variable_name, self.x, self.y,
                                       self.recon_grid, self.variable_units)

            # retrieve model data on identical grid, assign interpolators
            self.model_data = self.model_grid()
            self.create_2Dinterpolator(self.variable_name+'_model', self.x,
                                       self.y, self.model_data,
                                       self.variable_units)
            pdiff = (self.model_data - self.recon_grid)/self.model_data * 100.
            self.create_2Dinterpolator('PercentDiff', self.x, self.y,
                                       pdiff, '')
        elif self.run_option == 'flythrough':
            # fly through model data, accounting for shifts and wrapping of lon
            results = self.flythrough()  # full dictionary from flythrough
            # sort data onto a grid, assign interpolator
            self.recon_grid = self.grid_data(results)
            self.create_2Dinterpolator(self.variable_name, self.x, self.y,
                                       self.recon_grid, self.variable_units)
        else:
            raise AttributeError(f'Run option unknown: {self.run_option}.')
        print('Reconstruction program complete in ' +
              f'{ti.perf_counter()-t_start:.5f} s.')
        return

    # ---------------------- begin functions ----------------------------------

    def flythrough(self):
        '''fly trajectories through model data.'''

        # determine max value of longitude from coordinate type if spherical
        if self.coord_grid == 'sph' and self.coord_type == 'SPH':
            self.lon_max = 360.
        elif self.coord_grid == 'sph' and self.coord_type != 'SPH':
            self.lon_max = 180.

        # calculate offset trajectories and assemble into one array per coord
        fly_t, fly_c1, fly_c2, fly_c3 = [], [], [], []
        for i in range(len(self.time_offsets)):
            # calculate new trajectory
            new_time = self.sat_time + self.time_offsets[i]
            new_c1 = self.c1 + self.c1_offsets[i]
            new_c2 = self.c2 + self.c2_offsets[i]
            new_c3 = self.c3 + self.c3_offsets[i]

            # perform corrections and checks for lon and lat if in sph coords
            if self.coord_grid == 'sph':
                if self.c2_offsets[i] > 0 and (new_c2.min() < -90. or
                                               new_c2.max() > 90.):
                    print('Cannot correct for latitude shifts resulting in ' +
                          'latitudes greater than 90 degrees or less than ' +
                          '-90 degrees. Consider decreasing or removing the ' +
                          f'latitude offset {self.c2_offsets[i]}.')

                # check for and perform any longitude wrapping corrections
                while new_c1.max() > self.lon_max:
                    idx = np.where(new_c1 > self.lon_max)[0]
                    new_c1[idx] -= 360.

            # append arrays to list
            fly_t.append(new_time)
            fly_c1.append(new_c1)
            fly_c2.append(new_c2)
            fly_c3.append(new_c3)

        # change into 1D arrays
        fly_t = np.ravel(np.array(fly_t))
        fly_c1 = np.ravel(np.array(fly_c1))
        fly_c2 = np.ravel(np.array(fly_c2))
        fly_c3 = np.ravel(np.array(fly_c3))
        self.fly_t, self.fly_c1, self.fly_c2, self.fly_c3 = fly_t, fly_c1,\
            fly_c2, fly_c3
        ll = len(fly_t)

        # store average of entire trajectory if needed next or later or both
        if 'Avg' in self.recon_option:
            if self.t_avg == '':
                self.t_avg = np.mean(fly_t)
            if self.c1_avg == '':
                self.c1_avg = np.mean(fly_c1)
            if self.c2_avg == '':
                self.c2_avg = np.mean(fly_c2)
            if self.c3_avg == '':
                self.c3_avg = np.mean(fly_c3)

        # if AvgMod_AvgSlice is chosen, calculate avg of dims not reconstructed
        if self.recon_option[:6] == 'AvgMod':
            if 't' not in self.recon_dims:
                fly_t = np.repeat(self.t_avg, ll)
                print('Average time used for trajectories:', self.t_avg)
            if 'c1' not in self.recon_dims:
                fly_c1 = np.repeat(self.c1_avg, ll)
                print('Average c1 used for trajectories:', self.c1_avg)
            if 'c2' not in self.recon_dims:
                fly_c2 = np.repeat(self.c2_avg, ll)
                print('Average c2 used for trajectories:', self.c2_avg)
            if 'c3' not in self.recon_dims:
                fly_c3 = np.repeat(self.c3_avg, ll)
                print('Average c3 used for trajectories:', self.c3_avg)

        if 'OrbitSlice' in self.recon_option:
            # retrieve local times through SM coordinate conversion
            print('Filtering constellation trajectories on day/night option ' +
                  f'using GSE cartesian coordinates for {fly_t.size} ' +
                  'locations...')
            from kamodo_ccmc.flythrough.utils import ConvertCoord
            c1_GSE, c2_GSE, c3_GSE, units = ConvertCoord(
                fly_t, fly_c1, fly_c2, fly_c3, self.coord_type,
                self.coord_grid, 'GSE', 'car')

            # split into day or night based on option chosen
            # In GSE sph coordinates, day is lon>90 or lon<-90
            if 'OrbitSliceD' in self.recon_option:
                LT_idx = np.where(c1_GSE >= 0.)[0]  # x points to sun
            elif 'OrbitSliceN' in self.recon_option:
                LT_idx = np.where(c1_GSE <= 0.)[0]
            fly_t, fly_c1 = fly_t[LT_idx], fly_c1[LT_idx]
            fly_c2, fly_c3 = fly_c2[LT_idx], fly_c3[LT_idx]

        # decrease computation time by removing times not in files
        beg, end = MW.File_Times(self.model, self.file_dir, print_output=False)
        beg_ts, end_ts = beg.timestamp(), end.timestamp()
        tidx = np.where((fly_t >= beg_ts) & (fly_t <= end_ts))[0]
        if len(tidx) > 0:
            if len(tidx) < len(fly_t):
                diff = len(fly_t) - len(tidx)
                print(f'{diff} times not in model data. Removing.')
            fly_t, fly_c1 = fly_t[tidx], fly_c1[tidx]
            fly_c2, fly_c3 = fly_c2[tidx], fly_c3[tidx]

        # perform flythrough
        print('Performing satellite constellation flythrough for ' +
              f'{fly_t.size} locations...')
        results = SF.ModelFlythrough(
            self.model, self.file_dir, [self.variable_name], fly_t, fly_c1,
            fly_c2, fly_c3, self.coord_type+'-'+self.coord_grid,
            _print_units=False)

        # save correct time array for flythrough (all bad times removed)
        self.time_min = results['utc_time'].min()
        self.time_max = results['utc_time'].max()
        if self.t_avg > self.time_max or self.t_avg < self.time_min:
            self.t_avg = np.mean([self.time_min, self.time_max])
            print('Average time not in model data. Choosing new time in ' +
                  f'middle of model data time range: {self.t_avg}s or ' +
                  f'{ts_to_dt(self.t_avg)}')
        self.data_datetime = ts_to_dt(results['utc_time'][0])
        # self.results = results  # too much demand on memory
        return results

    def grid_data(self, results):
        '''Grid data along axes indicated. Handles all six options.'''
        t0 = ti.perf_counter()

        # collect variable data from satellite flythrough results
        variable_data = results[self.variable_name]

        # perform specialized gridding for sph coords (if lat or lon involved)
        if self.recon_dims == 'c1c2' and self.coord_grid == 'sph':  # 'LonLat'
            # initialize data
            x_data = results['c1']
            y_data = results['c2']

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create lon, lat coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int(360./self.dx) + 1
                self.ny = int(180./self.dy) + 1
                # center positions
                self.x = np.linspace(self.lon_max-360., self.lon_max, self.nx)
                self.y = np.linspace(-90., 90., self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Gridding data for LonLat...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        elif self.recon_dims == 'c2c3' and self.coord_grid == 'sph':
            # initialize data and height boundary
            x_data = results['c2']
            y_data = results['c3']

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create lat, height coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int(180./self.dx) + 1
                self.ny = int((y_data.max()-y_data.min())/self.dy) + 1
                self.x = np.linspace(-90., 90., self.nx)  # center positions
                self.y = np.linspace(y_data.min(), y_data.max(), self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Gridding data for LatHR...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        elif self.recon_dims == 'c1c3' and self.coord_grid == 'sph':
            # initialize data and height boundary
            x_data = results['c1']
            y_data = results['c3']

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create lon, height coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int(360./self.dx) + 1
                self.ny = int((y_data.max()-y_data.min())/self.dy) + 1
                # center positions
                self.x = np.linspace(self.lon_max-360., self.lon_max, self.nx)
                self.y = np.linspace(y_data.min(), y_data.max(), self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Gridding data for LonHR...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        elif self.recon_dims == 'tc2' and self.coord_grid == 'sph':
            # initialize data and time boundary
            xutc_data = results['utc_time']
            # convert to hours for plotting
            x_data = (xutc_data-xutc_data.min()) / 3600.
            y_data = results['c2']

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create time, lat coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int((x_data.max()-x_data.min())/self.dx) + 1
                self.ny = int(180./self.dy) + 1
                # center positions
                self.x = np.linspace(x_data.min(), x_data.max(), self.nx)
                self.y = np.linspace(-90., 90., self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Gridding data for TimeLat...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        elif self.recon_dims == 'tc1' and self.coord_grid == 'sph':
            # initialize data and time boundary
            xutc_data = results['utc_time']
            # convert to hours for plotting
            x_data = (xutc_data-xutc_data.min()) / 3600.
            y_data = results['c1']

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create time, lon coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int((x_data.max()-x_data.min())/self.dx) + 1
                self.ny = int(360./self.dy) + 1
                # center positions
                self.x = np.linspace(x_data.min(), x_data.max(), self.nx)
                self.y = np.linspace(self.lon_max-360., self.lon_max, self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Gridding data for TimeLon...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        elif 't' in self.recon_dims:  # 'tc3' sph and all car time options
            # initialize data and boundaries
            xutc_data = results['utc_time']
            # convert to hours for plotting
            x_data = (xutc_data-xutc_data.min()) / 3600.
            y_name = self.recon_dims[-2:]
            y_data = results[y_name]

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create time, height coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int((x_data.max()-x_data.min())/self.dx) + 1
                self.ny = int((y_data.max()-y_data.min())/self.dy) + 1
                # center positions
                self.x = np.linspace(x_data.min(), x_data.max(), self.nx)
                self.y = np.linspace(y_data.min(), y_data.max(), self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Gridding data for TimeHR...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        else:  # cartesian
            # initialize data and boundaries
            x_name = self.recon_dims[:-2]
            x_data = results[x_name]
            y_name = self.recon_dims[-2:]
            y_data = results[y_name]

            # define grid if not already defined
            if not (hasattr(self, 'x') and hasattr(self, 'y')):
                # create coordinate grid
                self.xmin, self.xmax = x_data.min(), x_data.max()
                self.ymin, self.ymax = y_data.min(), y_data.max()
                self.nx = int((x_data.max()-x_data.min())/self.dx) + 1
                self.ny = int((y_data.max()-y_data.min())/self.dy) + 1
                # center positions
                self.x = np.linspace(x_data.min(), x_data.max(), self.nx)
                self.y = np.linspace(y_data.min(), y_data.max(), self.ny)
            else:
                print('Using predefined grid.')
            variable = np.zeros((self.nx, self.ny))

            # sort data into grid
            print('Building cartesian reconstruction grid...', end="")
            for i in range(self.nx):
                for j in range(self.ny):
                    idx = np.where((abs(x_data-self.x[i]) < self.dx/2.) &
                                   (abs(y_data-self.y[j]) < self.dy/2.))[0]
                    if len(idx) > 0:
                        variable[i, j] = np.nanmean(variable_data[idx])
                    else:
                        variable[i, j] = np.nan

        print(f'done in {ti.perf_counter()-t0:.5f}s for {self.nx*self.ny} ' +
              'gridpoints.\n')
        return variable

    def model_grid(self):
        '''Retrieve model grid in method chosen.'''

        # section off code by analysis choice
        # UnMod_AvgSlice, AvgMod_AvgSlice
        if self.recon_option[-8:] == 'AvgSlice':
            # get slice from model data at avg value of each dim not used
            # determine avg value of dimensions not included in reconstruction
            if 't' not in self.recon_dims:
                new_time = self.t_avg
                print('Average t used for model:', new_time)
            if 'c1' not in self.recon_dims:
                new_c1 = self.c1_avg
                print('Average c1 used for model:', new_c1)
            if 'c2' not in self.recon_dims:
                new_c2 = self.c2_avg
                print('Average c2 used for model:', new_c2)
            if 'c3' not in self.recon_dims:
                new_c3 = self.c3_avg
                print('Average c3 used for model:', new_c3)

            # assign values for model grid according to recon_dims chosen
            # x coordinate first
            if self.recon_dims[:-2] == 't':  # convert hrs to UTC timestamps
                new_time = self.x * 3600. + self.time_min  # model_x
            elif self.recon_dims[:-2] == 'c1':
                new_c1 = self.x  # model_x
            elif self.recon_dims[:-2] == 'c2':
                new_c2 = self.x  # model_x
            elif self.recon_dims[:-2] == 'c3':
                new_c3 = self.x  # model_x
            # then y coordinate. the y coordinate is never time
            if self.recon_dims[-2:] == 'c1':
                new_c1 = self.y  # model_y
            elif self.recon_dims[-2:] == 'c2':
                new_c2 = self.y  # model_y
            elif self.recon_dims[-2:] == 'c3':
                new_c3 = self.y  # model_y

            # create coordinate arrays of same length for flythrough
            tt, xx, yy, zz = np.meshgrid(new_time, new_c1, new_c2, new_c3)
            tt, xx, yy, zz = np.ravel(tt), np.ravel(xx), np.ravel(yy),\
                np.ravel(zz)

            # perform flythrough
            print(f'Performing averaged grid flythrough for {tt.size} ' +
                  'locations...')
            t0 = ti.perf_counter()
            model_results = SF.ModelFlythrough(
                self.model, self.file_dir, [self.variable_name], tt, xx, yy,
                zz, self.coord_type+'-'+self.coord_grid, _print_units=False)
            # self.model_results = model_results  # store for sanity checks
            print(f'Grid flythrough completed in {ti.perf_counter()-t0:.5f} ' +
                  's.')

            # reshape data and return. The times should always be defined.
            if self.recon_dims == 'tc1':  # needed only for this combination
                model_data = np.reshape(model_results[self.variable_name],
                                        (self.y.size, self.x.size)).T
            else:
                model_data = np.reshape(model_results[self.variable_name],
                                        (self.x.size, self.y.size))
            return model_data

        # A model orbit slice is what an infinite number of satellites would
        # see if distributed along the orbit slice.
        # The longitude logic here compensates to keep equatorial crossing in
        # GSE at the same MLT as the original set.
        elif 'OrbitSlice' in self.recon_option:
            print('Performing orbit slicing of model data...')
            # determine approx length of an orbit
            lat_diff = np.diff(self.c2)  # find dlat
            idx = np.array([i for i in range(len(lat_diff)-1) if
                            lat_diff[i] > 0 and lat_diff[i+1] < 0],
                           dtype=int) + 2  # where orbits start
            T = []  # list to store orbital periods in for averaging later
            for i in range(len(idx)-1):  # ignore last orbit since partial
                singleorbit_time = self.sat_time[idx[i]: idx[i+1]]  # time vals
                T.append(singleorbit_time.max() - singleorbit_time.min())
            # find average of periods to nearest second
            orbital_period = int(np.mean(np.array(T)))
            print(f'Approximate orbital period is {orbital_period} s.')

            # initialize time offset grid and corresponding lon offsets
            tgrid = np.linspace(0, orbital_period,
                                int(orbital_period/self.dt) + 1)
            time_grid = np.repeat(tgrid, len(self.fly_t))
            lon_grid = -time_grid / 3600. * 15.  # -15 deg lon for each hour

            # apply grid of offsets to full constellation trajectory
            model_t = np.tile(self.fly_t, len(tgrid)) + time_grid
            model_c1 = np.tile(self.fly_c1, len(tgrid)) + lon_grid
            model_c2 = np.tile(self.fly_c2, len(tgrid))
            model_c3 = np.tile(self.fly_c3, len(tgrid))

            # retrieve local times through SM coordinate conversion
            print('Filtering model orbit slices on day/night option using ' +
                  f'GSE cartesian coordinates for {model_t.size} locations...')
            from kamodo_ccmc.flythrough.utils import ConvertCoord
            c1_MAG, c2_MAG, c3_MAG, units = ConvertCoord(
                model_t, model_c1, model_c2, model_c3, self.coord_type,
                self.coord_grid, 'GSE', 'car')

            # split into day or night based on option chosen
            # In GSE sph coordinates, day is lon>90 or lon<-90
            # x points to sun, x>=0 is dayside
            if 'OrbitSliceD' in self.recon_option:
                LT_idx = np.where(c1_MAG >= 0.)[0]
            elif 'OrbitSliceN' in self.recon_option:
                LT_idx = np.where(c1_MAG <= 0.)[0]
            model_t = model_t[LT_idx]
            model_c1 = model_c1[LT_idx]
            model_c2 = model_c2[LT_idx]
            model_c3 = model_c3[LT_idx]

            # perform flythrough
            print(f'{model_t.size} locations remaining for flythrough.')
            t0 = ti.perf_counter()
            model_results = SF.ModelFlythrough(
                self.model, self.file_dir, [self.variable_name], model_t,
                model_c1, model_c2, model_c3,
                self.coord_type+'-'+self.coord_grid, _print_units=False)

            # the orbit slices and day/night filtering result in orbit arrays
            # of different lengths, so can't just reshape. Have to regrid.
            model_data = self.grid_data(model_results)
            return model_data

        else:  # find average of slices in range of dimensions not used
            # assign values for model grid according to recon_dims option
            # chosen. *_AvgDSlice. x coordinate first.
            if self.recon_dims[:-2] == 't':  # convert hrs to UTC timestamps
                new_time = self.x*3600. + self.time_min  # model_x
            elif self.recon_dims[:-2] == 'c1':
                new_c1 = self.x  # model_x
            elif self.recon_dims[:-2] == 'c2':
                new_c2 = self.x  # model_x
            elif self.recon_dims[:-2] == 'c3':
                new_c3 = self.x  # model_x
            # then y coordinate. the y coordinate is never time
            if self.recon_dims[-2:] == 'c1':
                new_c1 = self.y  # model_y
            elif self.recon_dims[-2:] == 'c2':
                new_c2 = self.y  # model_y
            elif self.recon_dims[-2:] == 'c3':
                new_c3 = self.y  # model_y

            # apply d1 and d2 to other two dimensions to decrease computation
            # time. The t and c1 dimensions need to be switched in order for
            # the data to get reshaped correctly. It simply doesn't work the
            # normal way in 4D.
            if self.recon_dims == 'tc1':
                self.n1 = int(np.ceil((self.c2.max()-self.c2.min())/self.d1))
                self.n2 = int(np.ceil((self.c3.max()-self.c3.min())/self.d2))
                new_c2 = np.linspace(self.c2.min(), self.c2.max(), self.n1)
                new_c3 = np.linspace(self.c3.min(), self.c3.max(), self.n2)
                # flip t and c1
                final_shape = (self.ny, self.nx, self.n1, self.n2)
                mean_axes = (2, 3)
            elif self.recon_dims == 'tc2':
                self.n1 = int(np.ceil((self.c1.max()-self.c1.min())/self.d1))
                self.n2 = int(np.ceil((self.c3.max()-self.c3.min())/self.d2))
                new_c1 = np.linspace(self.c1.min(), self.c1.max(), self.n1)
                new_c3 = np.linspace(self.c3.min(), self.c3.max(), self.n2)
                # flip t and c1
                final_shape = (self.n1, self.nx, self.ny, self.n2)
                mean_axes = (1, 3)
            elif self.recon_dims == 'tc3':
                self.n1 = int(np.ceil((self.c1.max()-self.c1.min())/self.d1))
                self.n2 = int(np.ceil((self.c2.max()-self.c2.min())/self.d2))
                new_c1 = np.linspace(self.c1.min(), self.c1.max(), self.n1)
                new_c2 = np.linspace(self.c2.min(), self.c2.max(), self.n2)
                # flip t and c1
                final_shape = (self.n1, self.nx, self.n2, self.ny)
                mean_axes = (1, 2)
            elif self.recon_dims == 'c1c2':
                self.n1 = int(np.ceil((self.time_max-self.time_min)/self.d1))
                self.n2 = int(np.ceil((self.c3.max()-self.c3.min())/self.d2))
                new_time = np.linspace(self.time_min, self.time_max, self.n1)
                new_c3 = np.linspace(self.c3.min(), self.c3.max(), self.n2)
                # have to flip t and c1
                final_shape = (self.nx, self.n1, self.ny, self.n2)
                mean_axes = (0, 3)
            elif self.recon_dims == 'c1c3':
                self.n1 = int(np.ceil((self.time_max-self.time_min)/self.d1))
                self.n2 = int(np.ceil((self.c2.max()-self.c2.min())/self.d2))
                new_time = np.linspace(self.time_min, self.time_max, self.n1)
                new_c2 = np.linspace(self.c2.min(), self.c2.max(), self.n2)
                # flip t and c1
                final_shape = (self.nx, self.n1, self.n2, self.ny)
                mean_axes = (0, 2)
            elif self.recon_dims == 'c2c3':
                self.n1 = int(np.ceil((self.time_max-self.time_min)/self.d1))
                self.n2 = int(np.ceil((self.c1.max()-self.c1.min())/self.d2))
                new_time = np.linspace(self.time_min, self.time_max, self.n1)
                new_c1 = np.linspace(self.c1.min(), self.c1.max(), self.n2)
                # flip t and c1
                final_shape = (self.n2, self.n1, self.nx, self.ny)
                mean_axes = (0, 1)
            final_shape = list(final_shape)

            # create coord arrays of same length, perform flythrough
            tt, xx, yy, zz = np.meshgrid(new_time, new_c1, new_c2, new_c3,
                                         indexing='xy')
            print('Performing non-averaged grid flythrough for ' +
                  f'{np.ravel(tt).size} locations.')
            t0 = ti.perf_counter()
            model_results = SF.ModelFlythrough(
                self.model, self.file_dir, [self.variable_name], np.ravel(tt),
                np.ravel(xx), np.ravel(yy), np.ravel(zz),
                self.coord_type+'-'+self.coord_grid, _print_units=False)
            print(f'Grid flythrough completed in {ti.perf_counter()-t0:.5f} ' +
                  's.')

            # reshape, avg grids over non-reconstructed dimensions and return
            tmp = np.reshape(model_results[self.variable_name],
                             tuple(final_shape))
            model_array = np.transpose(tmp, axes=(1, 0, 2, 3))
            # ignore NaN values in avg
            model_data = np.nanmean(model_array, mean_axes)
            return model_data

    def xvec_dict(self):
        '''determine argument units based on dimension and coordinate choices
        '''

        xvec_dependencies = {}
        if 't' in self.recon_dims:
            xvec_dependencies['Time'] = 'hr'
        if 'c1' in self.recon_dims:
            if self.coord_grid == 'car':
                xvec_dependencies['X'] = 'R_E'
            elif self.coord_grid == 'sph':
                xvec_dependencies['Lon'] = 'deg'
        if 'c2' in self.recon_dims:
            if self.coord_grid == 'car':
                xvec_dependencies['Y'] = 'R_E'
            elif self.coord_grid == 'sph':
                xvec_dependencies['Lat'] = 'deg'
        if 'c3' in self.recon_dims:
            if self.coord_grid == 'car':
                xvec_dependencies['Z'] = 'R_E'
            elif self.coord_grid == 'sph' and self.coord_type == 'GDZ':
                xvec_dependencies['Alt'] = 'km'
            elif self.coord_grid == 'sph' and self.coord_type != 'GDZ':
                xvec_dependencies['Radius'] = 'R_E'
        return xvec_dependencies

    def create_2Dinterpolator(self, varname, x, y, data, units):
        '''Create gridded interpolator for reconstruction for any six dimension
        pairs.'''

        if x.size == y.size:
            data = data.T  # interpolator gets flipped somehow if square array
        xvec_dependencies = self.xvec_dict()

        # create gridded interpolator, use nearest-neighbor to avoid problems
        # with NaNs.
        rgi = RegularGridInterpolator((x, y), data, bounds_error=False,
                                      fill_value=None, method='nearest')

        '''#generic case
        @kamodofy(units = units, data = data)
        @gridify(x=x, y=y)
        def interpolator_grid(xvec):
            return rgi(xvec)
        '''
        # gridded interpolator inputs must be named according to variables
        # variables not allowed, so need text instead
        coord_list = list(xvec_dependencies.keys())
        if coord_list == ['Time', 'Lon']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, Lon=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Time', 'Lat']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, Lat=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Time', 'Alt']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, Alt=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Time', 'Radius']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, Radius=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Time', 'X']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, X=y)
            def interpolator_grid(xvec):
                """Interpolates 3d variable into a grid"""
                return rgi(xvec)

        elif coord_list == ['Time', 'Y']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, Y=y)
            def interpolator_grid(xvec):
                """Interpolates 3d variable into a grid"""
                return rgi(xvec)

        elif coord_list == ['Time', 'Z']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Time=x, Z=y)
            def interpolator_grid(xvec):
                """Interpolates 3d variable into a grid"""
                return rgi(xvec)

        elif coord_list == ['Lon', 'Lat']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Lon=x, Lat=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Lon', 'Alt']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Lon=x, Alt=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Lon', 'Radius']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Lon=x, Radius=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Lat', 'Alt']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Lat=x, Alt=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['Lat', 'Radius']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Lat=x, Radius=y)
            def interpolator_grid(thetavec):
                """Interpolates 3d variable into a grid"""
                return rgi(thetavec)

        elif coord_list == ['X', 'Y']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(X=x, Y=y)
            def interpolator_grid(xvec):
                """Interpolates 3d variable into a grid"""
                return rgi(xvec)

        elif coord_list == ['X', 'Z']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(X=x, Z=y)
            def interpolator_grid(xvec):
                """Interpolates 3d variable into a grid"""
                return rgi(xvec)

        elif coord_list == ['Y', 'Z']:
            @kamodofy(units=units, data=data, arg_units=xvec_dependencies)
            @gridify(Y=x, Z=y)
            def interpolator_grid(xvec):
                """Interpolates 3d variable into a grid"""
                return rgi(xvec)
        else:
            print(coord_list)
            raise AttributeError('Coordinate combination not supported.')

        # store interpolator, units, data, and xvec dictionary and return
        self[varname] = interpolator_grid
        self.variables[varname] = dict(units=units, data=data,
                                       arg_units=xvec_dependencies)
        self.variables[varname]['xvec'] = xvec_dependencies
        return
