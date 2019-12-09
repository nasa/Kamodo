
# 2019-12-09 13:52:30.292970: clock-in

# 2019-12-09 13:50:43.574951: clock-out
* working on cli

# 2019-12-09 10:05:45.828781: clock-in

# 2019-12-05 21:55:15.506661: clock-out
* looking at Y-combinator for recursive anonymous functions

# 2019-12-05 21:54:44.630995: clock-in: T-90m

# 2019-11-27 12:42:29.395096: clock-out
* fixing cdf arrays

# 2019-11-27 12:06:51.013991: clock-in

# 2019-11-22 09:57:46.571462: clock-out

* regex matching https://stackoverflow.com/questions/1687620/regex-match-everything-but-specific-pattern

# 2019-11-22 09:09:13.743347: clock-in

# 2019-11-21 12:09:12.103191: clock-out

# 2019-11-21 11:01:35.147547: clock-in

# 2019-11-21 10:27:47.314091: clock-out
* cdflib: switching to multiindex for all dependencies

# 2019-11-21 09:50:46.684607: clock-in

# 2019-11-20 16:27:17.501485: clock-out
* moved docs/notebooks/kameleon/kameleon_gateway.py into kamodo/readers
* pushed recent work
* created pull request

# 2019-11-20 15:30:59.817396: clock-in

# 2019-11-20 12:42:32.541944: clock-out

* on forwarding defaults: consider leveraging a function's .data during composition.
this would allow each downstream function to have automatic defaults!
* need to push these changes

# 2019-11-20 12:06:16.296442: clock-in

# 2019-11-20 11:58:22.508485: clock-out
* Kamodofied cdflib!

# 2019-11-20 09:35:49.895291: clock-in

# 2019-11-13 09:40:03.974566: clock-out

# 2019-11-13 08:58:52.993108: clock-in

# 2019-11-12 20:37:52.561965: clock-out

# 2019-11-12 19:17:10.399952: clock-in

# 2019-11-12 18:55:01.500196: clock-out
* command line

# 2019-11-12 18:37:46.018928: clock-in

# 2019-11-11 10:20:48.198321: clock-out
* set up command-line plotting

# 2019-11-11 09:35:53.331821: clock-in

# 2019-11-08 16:51:57.754553: clock-out
* begin work on command line interface
* trying facebook's hydra cli architecture
* `pip install hydra-core --upgrade`

# 2019-11-08 15:53:02.859151: clock-in

# 2019-10-18 16:47:23.366813: clock-out
* looking at inverse mapping
* Kamodofied pytiegcm
* need to finish gridifying inverse mapping
* consider adding indexing option to gridify

# 2019-10-18 11:58:47.333047: clock-in

# 2019-10-17 17:54:11.880758: clock-out
* kamodofying pyTIEGCM

# 2019-10-17 17:34:05.273266: clock-in

# 2019-10-17 15:54:57.284128: clock-out
* Kamodofied pyTIEGCM

# 2019-10-17 13:20:49.585114: clock-in

# 2019-09-30 12:34:55.407880: clock-out: added kamodofied kameleon example
* Kameleon-kamodo bridge development
* Added kamodofied kameleon object

# 2019-09-30 10:00:33.528255: clock-in

# 2019-09-11 13:32:56.778768: clock-out
* Created Fieldline Tutorial
* made griddify return in `xy` indexing to work with map-to-plane
* added scatter plot to available plot_types
* Finished kamodofying_models tutorial
* Need to check if plotting map-to-plane works for each axis

# 2019-09-11 08:19:33.545586: clock-in

# 2019-09-10 14:49:49.344817: clock-out
* fixed cone colors in plotting

# 2019-09-10 11:39:05.825143: clock-in

# 2019-09-10 10:43:05.346628: clock-out
* removing nans from solver output in favor of MultiIndex seeds

# 2019-09-10 09:55:34.659931: clock-in

# 2019-09-09 13:43:18.711955: clock-out
* Added pandas i/o for 3d line and vector plots

# 2019-09-09 11:07:35.809509: clock-in

# 2019-09-09 10:52:08.229479: clock-out
* Added pointlike decorator

# 2019-09-09 09:37:01.946111: clock-in

# 2019-09-06 15:26:26.921477: clock-out
* Looking into LFM wrapper
* https://wiki.ucar.edu/display/LTR/pyLTR seems to have been written for python 2

# 2019-09-06 15:20:58.638859: clock-in

# 2019-09-05 17:04:36.774012: clock-out
* solver decorator
* dipole field test
* stopping integration at boundary
* can choose resolution of knots

# 2019-09-05 11:09:34.195206: clock-in

# 2019-09-04 15:29:02.593470: clock-out
* got fieldline tracer working (ivp solver)

# 2019-09-04 14:47:48.601879: clock-in

# 2019-09-04 10:57:19.685814: clock-out

# 2019-09-04 10:53:23.666033: clock-in

# 2019-09-04 10:00:42.819821: clock-out
* looking at complex parameters for streamlines

# 2019-09-04 09:42:08.493375: clock-in

# 2019-09-03 18:33:35.737108: clock-out
* working on streamlines

# 2019-09-03 16:44:51.693606: clock-in

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

* Had issues installing locally with `pip install .` Need to use `pip install -e .` instead.

# 2019-08-14 10:37:19.291831: clock-in

# 2019-08-12 12:14:17.649049: clock-out
* Adding example notebook on kamodofication

# 2019-08-12 10:27:35.850505: clock-in

# 2019-07-24 16:44:02.434437: clock-out
* Moving code into new repo
* Creating new kamodo-user environment from clean miniconda

# 2019-07-24 14:33:41.015853: clock-in

