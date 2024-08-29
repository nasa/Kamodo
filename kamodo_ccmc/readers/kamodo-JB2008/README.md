
# Kamodo-JB2008

This is a kamodo wrapper for the JB2008 model.

We are building on the Destopy library.

To get started, clone this repo

<!-- #region -->
```sh
git clone git@github.com:EnsembleGovServices/kamodo-JB2008.git
cd kamodo-JB2008
```
<!-- #endregion -->

# Building with compose


Build with compose:

<!-- #region -->
```sh
docker compose build jb2008
```
<!-- #endregion -->

```python
from jb.jb08 import KJB08

import datetime
import numpy as np
from plotly.offline import init_notebook_mode
import plotly.graph_objs as go
```

```python
kj = KJB08()
kj
```

```python
alt = 400
lon = 0
lat = -90

my_date = datetime.datetime(2023, 4, 1, 1, 8, 0, 0)

kj.rho(lon, lat, alt, my_date)
```

```python
dens = []
lons = np.linspace(-180, 180, 300)
for lon in lons:
    dens.append(kj.rho(lon, 70, alt, my_date))
    
init_notebook_mode(connected=True)
go.Figure(go.Scatter(x=lons, y=dens))
```

```python

```
