# 2019-09-03 14:09:22.239226: clock-out: T-6m
* Deduct 6m
* Handled pandas output for plotting 3d vector fields
* Fixing integration symbol for solver
* Fixed griddable vector fields

# 2019-09-03 13:07:54.277159: clock-in

# 2019-08-30 10:49:09.145324: clock-out

* Issue when calling util.parse_latex:

> ImportError: LaTeX parsing requires the antlr4 python package, provided by pip (antlr4-python2-runtime or antlr4-python3-runtime) or conda (antlr-python-runtime)

* reconsider `antlr4` dependency

# 2019-08-30 10:47:31.579538: clock-in

# 2019-08-28 16:30:55.501447: clock-out: added gridify decorator
* Expanding kamodofication tutorial to support more variables
* Looking at Numpy's generalized universal function API - would formalize our mapping to plotting routines
* Add example of point interpolation for trajectory

# 2019-08-28 10:24:56.645799: clock-in

# 2019-08-16 16:27:22.801411: clock-out: cleaned up imports, added Kamodofication example, bug fixes
* testing install without psutil
* creating conda package
* cleaned up imports
* added click dependency
* updated pip

# 2019-08-16 14:25:59.427481: clock-in

# 2019-08-14 16:18:12.261596: clock-out
* finishing kamodofied models example

* Had issues installing locally with `pip install .` Resolved by following tip [here](https://github.com/pypa/pip/issues/5247#issuecomment-410910018): 
```console
pip install --upgrade --force-reinstall pip==9.0.3
pip install xxx --disable-pip-version-check
pip install --upgrade pip
```

# 2019-08-14 10:37:19.291831: clock-in

# 2019-08-12 12:14:17.649049: clock-out
* Adding example notebook on kamodofication

# 2019-08-12 10:27:35.850505: clock-in

# 2019-07-24 16:44:02.434437: clock-out
* Moving code into new repo
* Creating new kamodo-user environment from clean miniconda

# 2019-07-24 14:33:41.015853: clock-in

