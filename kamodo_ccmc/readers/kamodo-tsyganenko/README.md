
# Kamodo-Tsyganenko

This is a kamodo wrapper for the Tsyganenko models.

We are building on the [PyGeopack](https://github.com/mattkjames7/PyGeopack) library from Matt James

To get started, clone this repo

<!-- #region -->
```sh
git clone https://github.com/nasa/Kamodo.git
cd Kamodo/kamodo_ccmc/readers/kamodo-tsyganenko
```
<!-- #endregion -->

## Build `ensemble/tsyganenko`

build with compose:

<!-- #region -->
```sh
docker compose build tsyganenko
```
<!-- #endregion -->

This should build and install all the dependencies for the wrapper.


## Initialization

The first time you run the container, you will need to initialize the data needed by the Tsyganenko models

<!-- #region -->
```sh
docker compose initialize
```
<!-- #endregion -->

The compose file will automatically create volumes on your host machine and mount them into the above directories:

```yaml
volumes:
  - kp_data:/data/kp/
  - omni_data:/data/omni/
  - geopack_data:/data/geopack/
```



## Usage

```sh
from kamodo_geopack.tsyganenko import KTsyganenko, to_utc
```

```sh
kt = KTsyganenko(coord_out='GSM', coord_in='GSE')
kt
```

## Bfield

```sh
import pandas as pd
import numpy as np
```

To generate magnetic field vectors from components

```sh
t1 = pd.Timestamp.now() - pd.Timedelta('30d')
```

```sh
kt.Bvec(-30, 0, 0, t1)
```

```sh
kt.Bvec_n3(np.array([-30, 0, 0]), t1)
```

## Model Drivers

```sh
kt.V(t1)
```

```sh
kt.K_p(t1)
```

## Bfield trace


First we'll construct a set of seed points to pass to the functionalized tracer

```sh
import numpy as np

# Number of points
num_points = 30

# Radius of the circle
radius = 9.7

# Angles for each point
angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)

# Coordinates
x_coords = radius * np.cos(angles)
y_coords = radius * np.sin(angles)
z_coords = np.zeros(num_points)  # All z coordinates are 0
```

```sh
points = np.vstack((x_coords, y_coords, z_coords)).T
points.shape
```

Next choose a time for the model to trace

```sh
import pandas as pd
```

```sh
t1 = pd.Timestamp.now()
```

```sh
f = next(kt.F_B(points, t=t1-pd.Timedelta('100d')))
```

```sh
f
```

### Function oriented trace

To generate a collection of field line solutions from input seeds

```sh
fieldlines = kt.F_B(svec=points, t=t1-pd.Timedelta('100d'))
```

The above line should execute instantly because the result is a generator of functions.


You may iterate over each function, handling the traces on a case by case basis

```sh
f = next(fieldlines) # returns a function representing the field line through the next seed point
f
```

```sh
f().shape # calling this function with no arguments returns a fixed list of points 
```

The above points are linearly interpolated from the raw field line.


If you want to generate all field lines at once, construct a list.


```sh
fieldlines = list(fieldlines)
```

```sh
fieldlines
```

```sh
fieldlines[0]
```

### Object oriented tracing


If you prefer to work with the underlying field line objects directly, you can access the geopack api for the trace.


You will need to convert your time into date integer and float time components

```sh
from kamodo_geopack.tsyganenko import time_to_tsyg, gp
```

```sh
time = to_utc(t1)
date_int, ut = time_to_tsyg(time) # integer date, float time
```

```sh
date_int, ut
```

```sh
T = kt.trace(x=-30, y=0, z=0, date_int=date_int, ut=ut, coord_in='GSE')
```

access the positions as attributes

```sh
print(T.xgse.shape)
```

## Plotting traces



Plotting will iterate over all traces automatically.

```sh
from plotly.offline import init_notebook_mode
init_notebook_mode(connected=True)
```

```sh
fig = kt.plot(F_B=dict(svec=points, t=t1-pd.Timedelta('98d')))
# Set the range for x, y, z
fig.update_layout(
    scene=dict(
        xaxis=dict(range=[-30, 10]),  # Setting range for x-axis
        yaxis=dict(range=[-20, 20]),  # Setting range for y-axis
        zaxis=dict(range=[-20, 20]),  # Setting range for z-axis
#         aspectmode='manual',  # Allow manual setting for aspect ratio
#         aspectratio=dict(x=40, y=40, z=40)  # Aspect ratio corresponding to the ranges
    )
)
```

The magnetic field \(\vec{B}\) is defined as a vector field:
$$
\vec{B} = B_x(x, y, z) \vec{i} + B_y(x, y, z) \vec{j} + B_z(x, y, z) \vec{k}
$$

The unit vector of the magnetic field \(\hat{B}\) is given by:
$$
\hat{B} = \frac{\vec{B}}{|\vec{B}|} = \frac{B_x(x, y, z) \vec{i} + B_y(x, y, z) \vec{j} + B_z(x, y, z) \vec{k}}{\sqrt{B_x^2 + B_y^2 + B_z^2}}
$$

A curve \(\vec{r}(s)\), parametrized by arc length \(s\), traces the magnetic field lines if it satisfies the following differential equation:
$$
\frac{d\vec{r}}{ds} = \hat{B}(\vec{r}(s))
$$

This equation implies that the tangent to the curve at any point \(s\) is parallel and of the same magnitude (unit length) as the normalized magnetic field vector at that point. To solve for \(\vec{r}(s)\), integrate this equation:
$$
\vec{r}(s) = \vec{r}_0 + \int_{s_0}^{s} \hat{B}(\vec{r}(u)) \, du
$$

Here, \( \vec{r}_0 \) is the initial position vector, and \(s_0\) is the initial parameter value, typically representing the starting point of integration along the magnetic field line from a given seed point.



# Test commands

```sh
import PyGeopack as gp
```

```sh
x, y, z = -10, 0, 0
date = 20230411
ut = 5
```

<details>
<summary>Signature and Documentation for gp.ModelField Function</summary>
<pre>
Signature:
gp.ModelField(
    Xin,
    Yin,
    Zin,
    Date,
    ut,
    Model='T96',
    CoordIn='GSM',
    CoordOut='GSM',
    ReturnParams=False,
    **kwargs,
)
Docstring:
Calculates the model magnetic field at a given position or array of
positions in space.

Inputs
=======
Xin     : scalar or array containing the x positions(s).
Yin     : scalar or array containing the y positions(s).
Zin     : scalar or array containing the z positions(s).
Date    : Date - an integer in the format yyyymmdd.
ut      : Time in hours i.e. ut = h + m/60 + s/3600.
Model   : String to say which model to use out of the following:
        'T89'|'T96'|'T01'|'TS05' (see further below about models).
CoordIn : String denoting system of input coordinates out of:
        'GSE'|'GSM'|'SM'.
CoordOut        : String denoting output coordinate system out of:
        'GSE'|'GSM'|'SM'.

Keyword arguments
=================
iopt    : This keyword, if set, will override the iopt parameter
parmod  : set to a 10 element floating point array to override 
        the  default model params (see below for more info)
tilt    : Override for the actual dipole tilt angle (based on 
                the date and time) in radians
Vx              : X component solar wind velocity override.
Vy              : Y component solar wind velocity override.
Vz              : Z component solar wind velocity override.
Kp              : Sets the Kp index - essentially an override for iopt,
                where iopt = Kp + 1
Pdyn    : Override for parmod[0] - dynamic pressure
SymH    : Override for parmod[1] - SymH
By              : Override for parmod[2] - IMF By
Bz              : Override for parmod[3] - IMF Bz

Model Fields
============
'T89'   : The T89 model is only dependent upon the iopt parameter
                which is valid in the range 1 - 7, and is essentially
                equal to Kp + 1 (use iopt = 7 for Kp >= 6). No other 
                model uses the iopt parameter. Pdyn and Bz will be used to 
                check if we are inside the MP.
'T96'   : This model uses the first four parameters of parmod
                which are Pdyn, SymH, By and Bz, respectively. All other 
                elements of parmod are ignored.
'T01'   : This model uses the same 4 parameters as the T96 model,
                but also uses next two elements of the parmod array for
                the G1 and G2 parameters.
'TS05'  : This model uses all ten parmod elements - the first 4
                are the same as T96, the next 6 are the W1-W6 parameters
                which **apparently** aren't too important (you could
                probably set these to 0). The ones calculated by this 
                module are probably erroneous - they are calculated by
                using Tsyganenko's own Fortran code, but they still 
                don't match the ones provided on his website! So maybe
                you *should* set these to 0.
                
NOTE 1 : Vx, Vy and Vz are used to convert from GSM to GSW coords.
NOTE 2 : All parameters here will be automatically filled in 
        based on the OMNI data - any gaps could be replaced with 
        custom parameters.
NOTE 3 : When using custom parameters, all other parameters will
        still be calculated automatically, unless they are also
        customized.
</pre>
</details>


```sh
Bx,By,Bz = gp.ModelField(x, y, z, date, ut, Model='TS05',CoordIn='GSM',CoordOut='GSM')
```

<!-- #region -->
output:

```sh
[4.95974694] [2.185299] [13.73547959]

```
<!-- #endregion -->

```sh
T = gp.TraceField(x,y,z,date,ut)
```

```sh
from plotly.offline import init_notebook_mode
```

```sh
import plotly.graph_objs as go
```

```sh
init_notebook_mode(connected=False)
```

```sh
import numpy as np

# Number of points
num_points = 30

# Radius of the circle
radius = 9.7

# Angles for each point
angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)

# Coordinates
x_coords = radius * np.cos(angles)
y_coords = radius * np.sin(angles)
z_coords = np.zeros(num_points)  # All z coordinates are 0

x_, y_, z_ = np.array([]), np.array([]), np.array([])

# Print the coordinates
for x, y, z in zip(x_coords, y_coords, z_coords):
    T = gp.TraceField(x,y,z,date,ut)
    s = T.s[~np.isnan(T.s)]
    print(min(s), max(s))
    x_pnts = T.xgsm[0][~np.isnan(T.xgsm[0])]
    y_pnts = T.ygsm[0][~np.isnan(T.ygsm[0])]
    z_pnts = T.zgsm[0][~np.isnan(T.zgsm[0])]
    x_ = np.concatenate([x_, x_pnts, [np.nan]])
    y_ = np.concatenate([y_, y_pnts, [np.nan]])
    z_ = np.concatenate([z_, z_pnts, [np.nan]])
    
go.Figure([go.Scatter3d(x=x_, y=y_, z=z_, mode='lines')])
```

```sh
date
```

```sh
ut
```

```sh

```

```sh

```

```sh

```

```sh

```
