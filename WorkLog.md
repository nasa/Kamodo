
# 2020-11-17 16:06:38.651579: clock-in

# 2020-11-17 13:46:18.035319: clock-out

* passing functional unit test

# 2020-11-17 11:17:10.446484: clock-in

# 2020-11-16 18:34:07.584293: clock-out

* implementing unit functions

# 2020-11-16 16:59:30.697802: clock-in

# 2020-11-16 15:01:46.527485: clock-out


# 2020-11-16 13:11:21.254349: clock-in

# 2020-11-16 13:05:29.126719: clock-out: T-1h 


# 2020-11-16 11:12:27.592786: clock-in

# 2020-11-15 17:03:30.218928: clock-out


# 2020-11-15 15:17:08.285800: clock-in

# 2020-11-15 13:09:56.154786: clock-out


# 2020-11-15 11:52:34.355935: clock-in

# 2020-11-13 23:16:23.702830: clock-out


# 2020-11-13 22:37:29.395330: clock-in

# 2020-11-13 18:50:51.601791: clock-out


# 2020-11-13 17:01:21.677780: clock-in

# 2020-11-13 11:59:07.207912: clock-out

* overhauling units

# 2020-11-13 11:23:14.433529: clock-in

# 2020-11-12 15:25:13.259993: clock-out


# 2020-11-12 13:53:04.513345: clock-in

# 2020-11-12 11:40:29.694990: clock-out


# 2020-11-12 11:22:56.545479: clock-in

# 2020-11-11 17:42:08.055155: clock-out

* updated convert_to to raise errors

# 2020-11-11 16:46:46.677731: clock-in

# 2020-11-11 16:29:02.695711: clock-out


# 2020-11-11 16:10:38.700105: clock-in

# 2020-11-09 11:08:06.200728: clock-out

* fixing unit bugs

# 2020-11-09 10:00:15.420963: clock-in

* fixing unit conversion bugs

# 2020-11-03 11:17:20.736909: clock-out

* added dynamic function evaluation

# 2020-11-03 09:18:36.269712: clock-in

* radians
* why can't user access constants?

# 2020-10-27 23:26:44.585965: clock-out

* developing evaluate endpoint


`http://127.0.0.1:5000/api/mymodel/evaluate?variable=%27g=(f%2B1)**.5%27&x=[3,4,5]`

plus sign: `%2B`

# 2020-10-27 21:16:42.216480: clock-in

# 2020-10-27 15:30:20.449578: clock-out

* merging tests from Dhruv

# 2020-10-27 13:31:24.216328: clock-in

* code cleanup

# 2020-10-19 18:22:01.985135: clock-out

* all util.py unit tests pass
* spacing

* trying to fix collections warning

```bash
DeprecationWarning: Using or importing the ABCs from 'collections' instead of from 'collections.abc' is deprecated since Python 3.3,and in 3.9 it will stop working
```

# 2020-10-19 17:32:57.246898: clock-in

# 2020-10-14 23:22:57.047235: clock-out

* commenting util.py

# 2020-10-14 21:43:10.943127: clock-in

# 2020-10-14 12:31:50.544849: clock-out

* going through util.py

# 2020-10-14 11:29:55.853990: clock-in

# 2020-10-13 15:52:04.229990: clock-out

* cleaning up unit conversion code

# 2020-10-13 12:10:24.428950: clock-in

# 2020-09-02 16:03:51.164640: clock-out

* how to keep NASA readers and core from conflicting:
	- currently these are separate files so merges should be straight-forward
	- readers are all subclasses of Kamodo, so breaking changes should only be downstream
* modifications to core should be made with hourly to comply with NASA license
* could also rewrite kamodo core as a separate package with its own repo

# 2020-09-02 14:59:41.540200: clock-in

# 2020-08-12 12:54:08.590956: clock-out

* dev meeting

# 2020-08-12 11:59:10.698783: clock-in

# 2020-08-05 12:37:04.934436: clock-out

* dev meeting

# 2020-08-05 12:01:11.230431: clock-in

# 2020-07-29 12:58:25.948624: clock-out

* dev meeting
* modularity
	- make sure readers are independent
	- allow matplotlib without installing plotly

# 2020-07-29 12:02:08.702440: clock-in

# 2020-07-15 13:00:19.299408: clock-out

* dev meeting
* send blurb on cli to darren for gem

# 2020-07-15 12:31:56.828729: clock-in: T-15m 

* developer meeting
# 2020-07-08 12:31:42.478226: clock-out

* developer meeting
* attend GEM/AGU?

# 2020-07-08 11:59:12.669083: clock-in

# 2020-07-08 11:52:18.410732: clock-out


# 2020-07-08 11:50:08.658452: clock-in

# 2020-07-01 12:54:01.716646: clock-out

* developer meeting

# 2020-07-01 12:04:30.950750: clock-in

# 2020-06-24 12:54:16.183808: clock-out

* developer meeting

# 2020-06-24 12:16:33.594954: clock-in

# 2020-06-17 12:13:35.118654: clock-out


# 2020-06-17 12:12:50.788007: clock-in

# 2020-06-10 12:46:57.416636: clock-out

* developer meeting

# 2020-06-10 12:16:22.661881: clock-in

# 2020-06-03 13:01:51.004408: clock-out

* kamodo team meeting

# 2020-06-03 13:01:28.636263: clock-in: T-50m 

* forwarding plot args to plot funcs, works for quiver plots
* found forge tool to generate custom function signatures!

```console
pip install python-forge
```
```python
@forge.sign(
    forge.arg('x'),
    forge.arg('y'),
    forge.arg('opt', default=None),
)
def func(*args, **kwargs):
    # signature becomes: func(x, y, opt=None)
    return (args, kwargs)

assert func(1, 2) == ((), {'x': 1, 'y': 2, 'opt': None})
```

# 2020-05-27 14:38:38.650053: clock-out


# 2020-05-27 14:35:57.387894: clock-in

# 2020-05-27 12:57:54.286201: clock-out

* developer meeting

# 2020-05-27 11:38:16.287492: clock-in

# 2020-05-23 20:54:58.496791: clock-out

* fixed ordering bug in generator function evaluation

# 2020-05-23 20:01:46.123612: clock-in

* cleaning up skew contour carpet plots
* squeezing gridify output, added rvert lvert

# 2020-05-13 15:02:23.437971: clock-out

* got basic api working
* developer meeting
* making `get_defaults` return None for args without defaults
* argument parsing https://flask-restful.readthedocs.io/en/latest/reqparse.html
* got api to return model description 

* to specify curl `-d` options as url query, use `-G` flag

So this:

	curl http://127.0.0.1:8050/api/mymodel -d greeting=goodbye -G

is equivalent to

	http://127.0.0.1:8050/api/mymodel?greeting=goodbye

# 2020-05-13 09:06:37.493882: clock-in

# 2020-04-29 13:42:44.161370: clock-out

* flask server integration, api test working
* PYHC meeting

# 2020-04-29 12:14:17.372197: clock-in

# 2020-04-22 14:16:35.278555: clock-out: T-1h16m 

* More on flask integration from [plotly](https://dash.plotly.com/integrating-dash)


# 2020-04-22 11:56:00.335536: clock-in

# 2020-04-15 13:16:43.205737: clock-out: T-15m 


# 2020-04-15 11:57:50.421409: clock-in

# 2020-04-08 13:01:50.501409: clock-out

* meeting with developers

# 2020-04-08 12:18:16.752330: clock-in

# 2020-04-01 13:52:05.021374: clock-out

* read through [flask app tutorial](https://github.com/toddbirchard/plotlydash-flask-tutorial)
* ran with `python wsgi.py `

## Kamodo meeting

pandas time interpolation:

```python
def time_interpolation(df, t):
	# combine original time index with new times
	df_ = df.reindex(df.index.union(t))
	# apply pandas' interpolation method to fill in new data
	df_interpolated = df_.interpolate(method='time')
	# return only values at selected times
	result = df_interpolated.reindex(t)
	return result
```

plans for May:
* gui working for summer schools
* include field line tracing
* include services
* choose a name for kamodo network

# 2020-04-01 12:08:40.151618: clock-in

# 2020-03-25 14:02:55.608031: clock-out

## kamodo meeting
* set deadlines for gui, april->may
* need a function that converts rho(x,y,z) -> rho(xvec)
* [dash-flask app tutorial](https://hackersandslackers.com/plotly-dash-with-flask/)
* tried running flask_restful with Dash, got `'Dash' object has no attribute 'handle_exception'`

# 2020-03-25 12:30:23.986927: clock-in

# 2020-03-25 12:20:20.625681: clock-out


# 2020-03-25 12:04:20.154862: clock-in

# 2020-03-18 21:37:24.225992: clock-out: T-5h 


# 2020-03-18 16:26:22.789889: clock-in

# 2020-03-18 13:03:13.421786: clock-out

* looking at flask rest api https://flask-restful.readthedocs.io/en/latest/quickstart.html#full-example

## developer meeting
* sscweb is positions also extrapolated into the future
* cdaweb is positions and data

# 2020-03-18 11:31:13.223644: clock-in

# 2020-03-11 20:25:31.313266: clock-out

* developer meeting

# 2020-03-11 20:25:18.435898: clock-in: T-45m 

# 2020-03-05 00:02:00.846556: clock-out

* prototyping mas

# 2020-03-04 22:01:13.683543: clock-in: T-30m 

# 2020-03-04 13:12:38.209902: clock-out

## developer meeting
* put website url in github description
* pyHC standards grading

https://github.com/heliophysicsPy/heliophysicsPy.github.io/blob/master/_pyhc_projects/pyhc_project_grading_guidelines.md

# 2020-03-04 11:52:51.855025: clock-in

# 2020-02-25 15:49:05.616947: clock-out

* addressing Liang's suggestions
* fixed deprecation warning from sympy>=1.5

# 2020-02-25 14:37:23.285614: clock-in

* check out partial functions
# 2020-02-19 13:03:32.424355: clock-out

* developer meeting

# 2020-02-19 12:22:27.823834: clock-in

# 2020-02-16 17:09:14.266657: clock-out: T-56h 

* submitted several feature issues

# 2020-02-14 09:56:13.530302: clock-in

# 2020-02-13 10:39:58.363603: clock-out

* iSWAT-COSPAR sessions
* end-to-end solutions
* community involvement
* visibility
* understanding interpolators

# 2020-02-13 09:09:49.084779: clock-in

# 2020-02-13 09:09:44.638954: clock-out: T-20h 


# 2020-02-12 09:37:08.349586: clock-in: T-30m 

# 2020-02-10 12:05:07.684918: clock-out

* fixed cli bug that prevented multiple plots from being saved

# 2020-02-10 12:03:49.344691: clock-in: T-1h 

# 2020-02-09 22:22:51.749709: clock-out

* gui and cli release
* fixed continuous reload bug

# 2020-02-09 20:13:13.525126: clock-in

# 2020-02-09 11:07:16.234918: clock-out


# 2020-02-09 10:30:57.198701: clock-in

# 2020-02-08 19:06:04.754325: clock-out

* gui improvements

# 2020-02-08 17:41:35.076690: clock-in

# 2020-02-08 15:54:29.107857: clock-out

* got reload config to work\!

# 2020-02-08 15:17:25.461178: clock-in

# 2020-02-08 14:26:00.479448: clock-out

* tried getting fully dynamic dash callbacks to work

# 2020-02-08 12:40:38.991214: clock-in

# 2020-02-07 19:58:47.600209: clock-out

* got interactiv configuration

# 2020-02-07 18:28:09.035029: clock-in

# 2020-02-05 18:28:53.256342: clock-out


# 2020-02-05 16:56:46.743551: clock-in

# 2020-02-05 16:50:08.675833: clock-out

* testing stateful storage update

# 2020-02-05 16:22:00.871936: clock-in

# 2020-02-05 14:14:54.356827: clock-out

* got clientside subplots to render
* meeting for iSWAT-COSPAR
* pushed code into NASA master
* need to update pypi version
* think about kamodo api
* docker container!

# 2020-02-05 11:23:08.377881: clock-in

# 2020-02-04 14:59:50.123294: clock-out

* got gui to load separate models and parameters

# 2020-02-04 13:22:32.543080: clock-in

# 2020-02-04 11:59:54.684191: clock-out

* reading up on dcc.store and clientside_callback

# 2020-02-04 10:39:40.567881: clock-in

# 2020-02-03 16:29:29.417205: clock-out

* added parameter checkboxes

# 2020-02-03 15:49:51.793947: clock-in

# 2020-02-03 12:36:08.959116: clock-out

* rendering equations through katex

# 2020-02-03 10:22:19.591903: clock-in

# 2020-01-31 19:35:30.915089: clock-out

* got range slider to work

# 2020-01-31 18:37:19.539544: clock-in

# 2020-01-31 17:56:33.727653: clock-out


# 2020-01-31 17:34:13.290486: clock-in

# 2020-01-31 17:19:16.956235: clock-out


# 2020-01-31 17:02:15.856951: clock-in

# 2020-01-31 16:52:07.773236: clock-out


# 2020-01-31 16:06:35.834953: clock-in

# 2020-01-31 15:11:15.638982: clock-out


# 2020-01-31 14:49:48.788349: clock-in

# 2020-01-31 13:27:55.782459: clock-out

* packaging pysatKamodo

# 2020-01-31 12:37:50.196753: clock-in

# 2020-01-31 11:19:55.778624: clock-out

* gui generating callbacks after layout is set

# 2020-01-31 10:35:50.687658: clock-in

# 2020-01-29 13:18:25.114169: clock-out

* meeting with Kamodo team for iSWAT-COSPAR prep

# 2020-01-29 11:59:08.216757: clock-in

* merging hourly.yaml
# 2020-01-29 11:19:29.007340: clock-out

* added link to NASA CCMC site
* cleaned up documentation site
* cleaning up docs

# 2020-01-29 09:50:29.055214: clock-in: T-20m 

# 2020-01-28 13:39:50.454679: clock-out

* answering support email

# 2020-01-28 13:30:23.612051: clock-in

# 2020-01-22 15:09:42.523706: clock-out


# 2020-01-22 15:09:36.423725: clock-in: T-40m 

* adding separate tabs for each model
# 2020-01-22 14:19:44.506406: clock-out

* working on gui layout

# 2020-01-22 12:13:09.639567: clock-in

# 2020-01-21 17:53:12.595685: clock-out

* added dynamic line plots to gui

# 2020-01-21 15:28:20.175086: clock-in

# 2020-01-21 13:37:38.575443: clock-out

* added composition to cli
* made compose use keyword arguments
* added compose function for multiple kamodo objects

# 2020-01-21 11:10:26.570097: clock-in

# 2020-01-20 15:23:39.237942: clock-out

* adding multiple models to cli
* setting up kamodo.yaml overrides
* Answered Lutz's issue

# 2020-01-20 11:21:44.788233: clock-in

# 2020-01-15 13:11:48.588900: clock-out

* testing hourly commit
* addressing issue from Lutz
* send agenda for meeting

# 2020-01-15 11:45:55.727665: clock-in

# 2019-12-20 13:50:20.467183: clock-out
* working on cli
* added ability to call kamodo from any work directory containing config.yaml
* added embedable plots: `model.plot_conf.output_type=div` saves the div to the file
* have markdown-include point to the div, so you can embed in your own docs!

# 2019-12-20 11:18:35.468801: clock-in

# 2019-12-19 12:17:51.613183: clock-out
* adding config file description to docs
* Looking at `Config search path`. 
* hydra plugin configs are found automatically once they are installed
* how do we use the cwd to configure? https://github.com/facebookresearch/hydra/issues/274
* docusaurus looks interesting

# 2019-12-19 11:26:13.027723: clock-in

# 2019-12-09 13:53:18.142633: clock-out
* forgot to include cli notebook

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

