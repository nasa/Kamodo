
### 2021-09-21 12:48:14.222004: clock-out

* comparing with kamodo-core, both have psi/ensemble branch
* the goal is to remove readers from kamodo-core and remove kamodo from ccmc. Need to make sure we're not losing any work

```sh
git remote add kamodo_core git@github.com:ensemblegov/kamodo-core.git
git fetch kamodo_core master:core_master # fetches master from kamodo_core and names it core_master branch
git diff core_master # compare master to core_master
```

* looking at packaging namespace packages https://packaging.python.org/guides/packaging-namespace-packages/
* a good blog post on py3.3+ namespace packaging https://newbedev.com/is-init-py-not-required-for-packages-in-python-3-3

### 2021-09-21 10:43:18.710093: clock-in: T-10m 

### 2021-09-20 17:41:43.068643: clock-out


### 2021-09-20 16:54:31.143729: clock-in

### 2021-09-20 11:50:02.590611: clock-out

* about section, created worklogs

### 2021-09-20 11:40:58.884429: clock-in

### 2021-06-21 16:18:17.807491: clock-out

* installing kamodo in editable mode
* building from source

### 2021-06-23 11:56:51.629840: clock-out


### 2021-06-23 11:01:39.298459: clock-in

### 2021-06-22 17:47:05.431045: clock-out

* moving animations into plot method

### 2021-06-22 16:30:22.893968: clock-in

### 2021-06-21 16:18:17.807491: clock-out

* using squeeze flag
* fixed test affected by size_threshold

### 2021-06-21 15:46:45.787242: clock-in

### 2021-06-21 13:55:39.710580: clock-out

* sped up animations

### 2021-06-21 13:13:21.460428: clock-in

### 2021-06-21 11:59:53.993441: clock-out

* adding animations - full_figure_for_development takes a long time to run. need a workaround

### 2021-06-21 11:35:38.562849: clock-in

### 2021-06-18 15:36:29.661450: clock-out

* got 3d animation working

### 2021-06-18 14:52:28.603029: clock-in

### 2021-06-18 11:22:00.529644: clock-out

* animating 2d parametric plots

### 2021-06-18 10:41:46.096653: clock-in

### 2021-06-17 16:35:09.484743: clock-out

* autoranged frames for 2d plots
* `pip install kaleido` needed to get frame layouts

### 2021-06-17 15:12:28.301751: clock-in

### 2021-06-17 12:34:51.152273: clock-out

* got animations working

### 2021-06-17 12:07:27.656044: clock-in

### 2021-06-16 16:09:33.197304: clock-out

* working on animations

### 2021-06-16 14:45:42.337094: clock-in

### 2021-06-16 12:16:27.943176: clock-out

* partial generators

### 2021-06-16 11:19:26.928668: clock-in

### 2021-06-15 19:42:28.601429: clock-out

* partial decorator passes functionality tests

### 2021-06-15 18:20:52.771686: clock-in

### 2021-06-15 14:39:28.250605: clock-out

* working on partial decorator

### 2021-06-15 13:19:46.744384: clock-in

### 2021-06-15 12:45:11.749450: clock-out

* problems with functools.partial decorator
Here is a problem with our curry operator is that we have to unpack all the arguments in order to evaluate

```python
@curry
def f(x=1,y=2,z=3):
    return x+y+z
```

`f(1)` returns a function of `y` instead of `y,z`, but that's how `currying` is [supposed to work](https://en.wikipedia.org/wiki/Currying).

To return a function of two variables, we are back to partials. It seems that these are separate decorators.


### 2021-06-15 12:43:53.208925: clock-in: T-41m 


### 2021-06-14 20:23:51.337807: clock-out


### currying decorator - kwargs

What we want is a decorator that returns a stateless function with a new signature like this:

```python
@curry
def f(x, y, z):
    return x + y + z

g = f(2)(3)

assert g(1) == 1 + 2 + 3
```



Normal python functions have `args` and `kwargs`:

* `args` are required
* `kwargs` are defaults

So one way to achieve the above behavior is to convert `args` into `kwargs`. `g(1)` would be equivalent to:

```python
g = lambda z, x=2, y=3: f(x, y, z)
g(1)
```

The problem is - what do we do with the original function defaults.


```python
@curry
def f(z=3, y=2, x=1):
    return x + y + z

f(1) == g(1)
```
Now currying has no effect. Recall why we are currying in the first place: we want to fix the values of a function so that we can evaluate it over a subset of the original arguments. If we already have defaults, do we really need to curry? We could have our curry operator eliminate the defaults? Or we could have the defaults pass through:

```python
g = curry(f)
assert g()()() == f()
```

This way we aren't losing any information from the original function, but are still allowing it to be called with single arguments.

```python
assert g()(1)() == f(y=1)
```

### 2021-06-14 19:31:40.201422: clock-in

### 2021-06-14 14:02:26.824598: clock-out


### 2021-06-14 13:46:07.459008: clock-in

### 2021-06-11 17:43:46.767514: clock-out

* trying out decorator curry

### 2021-06-11 16:05:44.670016: clock-in

### 2021-06-09 18:59:32.741722: clock-out


### 2021-06-09 18:45:08.502957: clock-in

### 2021-06-09 18:44:17.665337: clock-out: T-2h 

* added currying decorator
* [currying](https://www.python-course.eu/currying_in_python.php) in python.

### 2021-06-09 13:42:15.184948: clock-in: T-15m 

### 2021-06-09 10:02:25.693985: clock-out


### 2021-06-09 09:40:25.698556: clock-in

### 2021-06-08 19:33:26.008122: clock-out: T-10m 


### 2021-06-08 19:07:36.489347: clock-in

### 2021-06-08 18:07:48.196354: clock-out


### 2021-06-08 17:58:03.342088: clock-in

### 2021-06-08 17:33:41.236397: clock-out


### 2021-06-08 16:47:15.520867: clock-in

### 2021-06-07 19:18:28.125536: clock-out

* geneartor input type

### 2021-06-07 18:31:58.228494: clock-in

### 2021-06-02 14:52:38.570656: clock-out

* fixed bugs in plotting

### 2021-06-02 13:50:08.304688: clock-in

### 2021-06-02 10:00:40.767448: clock-out

* fixing plot bugs
* there is bug in equation rendering where Kamodo.to_latex returns the function's expression rather than the symbol used in registration.

### 2021-06-02 08:21:21.610958: clock-in

### 2021-06-01 19:41:20.602808: clock-out


### 2021-06-01 17:11:13.591391: clock-in

### 2021-06-01 13:52:37.320039: clock-out

* fixing plot key bugs

### 2021-06-01 12:48:58.347646: clock-in

### 2021-05-27 13:10:05.047450: clock-out

* looking at generator arguments

### 2021-05-27 12:40:51.097172: clock-in

### 2021-05-26 12:25:41.791573: clock-out

* issues with slice generation returning empty plot

### 2021-05-26 12:21:26.381522: clock-in

### 2021-05-25 19:47:12.438128: clock-out

* fixing pd.datetime deprecration warning

### 2021-05-25 19:39:05.600086: clock-in: T-20m 

### 2021-05-25 17:39:20.981299: clock-out

* made contour time slider

### 2021-05-25 17:31:43.349647: clock-in: T-30m 

### 2021-05-25 16:42:14.594993: clock-out


### 2021-05-25 15:48:09.688171: clock-in

### 2021-05-25 12:52:45.366801: clock-out

* flattening arg shapes in preparation for 4d

### 2021-05-25 12:09:08.836291: clock-in

### 2021-05-24 18:02:15.095789: clock-out

* refactored and added squeeze kwarg to gridify
* refactoring `@gridify` to use forge
* something to keep in mind when using forge:

```
.. warning::

    When supplying previously-created parameters to :func:`~forge.sign`,
    those parameters will be ordered by their creation order.

    This is because Python implementations prior to ``3.7`` don't
    guarantee the ordering of keyword-arguments.

    Therefore, it is recommended that when supplying pre-created
    parameters to :func:`~forge.sign`, you supply them as positional
    arguments:


    .. testcode::

        import forge

        param_b = forge.arg('b')
        param_a = forge.arg('a')

        @forge.sign(a=param_a, b=param_b)
        def func1(**kwargs):
            pass

        @forge.sign(param_a, param_b)
        def func2(**kwargs):
            pass

        assert forge.repr_callable(func1) == 'func1(b, a)'
        assert forge.repr_callable(func2) == 'func2(a, b)'

```


### 2021-05-24 16:55:07.640148: clock-in

* differential equations could be written as function decorators applied to boundary conditions.
* boundary conditions are functions over the simulation domain boundary with nans everywhere else.
* solutions are returned as functions over the independent variables. 

### 2021-05-21 16:47:38.435955: clock-out

* cleaned up plot key generation

### 2021-05-21 15:24:12.011043: clock-in

### 2021-05-21 13:07:05.626110: clock-out

* simplifying plot key generation

### 2021-05-21 11:12:40.628676: clock-in

look at jupyter book for publication https://jupyterbook.org/intro.html

### 2021-03-31 14:02:27.533428: clock-out

* registered image plot type

### 2021-03-31 13:02:51.467590: clock-in

* made kamodo unit system the default for conversion

### 2021-03-30 18:45:14.238802: clock-out

* working on functional images
* added angular frequency units

### 2021-03-30 16:16:32.599994: clock-in

### 2021-03-30 13:39:42.776447: clock-out


### 2021-03-30 12:21:01.381129: clock-in

### 2021-03-30 12:01:37.303279: clock-out


### 2021-03-30 11:11:21.672511: clock-in

### 2021-03-29 18:59:34.366365: clock-out

* adding plasmapy kamodofication bug test
* fixed bug in parse_expr locals
* fixed bug with Newton symbol clash

### 2021-03-29 17:01:31.209600: clock-in

### 2021-03-29 16:03:16.457691: clock-out

* getting function has no attribute subs with N - need to check against `_clash` list

### 2021-03-29 15:47:37.025035: clock-in

### 2021-03-29 12:04:33.031396: clock-out

* looking at images

### 2021-03-29 11:16:58.877386: clock-in

### 2021-03-27 14:11:10.779236: clock-out

* fixed example signature
* fixed workflow
* fixed bug in get_dimensions preventing certain unit conversions

### 2021-03-27 13:06:31.717147: clock-in

### 2021-03-27 12:52:56.983449: clock-out


### 2021-03-27 12:24:11.426813: clock-in

### 2021-03-27 11:58:55.516981: clock-out

* test of pascals failing
* fixed to_latex rendering

### 2021-03-27 10:42:18.972386: clock-in

### 2021-03-27 10:34:30.659703: clock-out: T-10m 

* cleariving cells
* fixed to-html
* adding notebooks list
* adding Visualization notebook to workflow

### 2021-03-27 09:32:23.833239: clock-in

### 2021-03-24 22:34:04.374828: clock-out

* fixed multi argument unit composition

### 2021-03-24 22:23:10.086953: clock-in

### 2021-03-24 21:42:25.386116: clock-out

* unify expr args not in same order as free_symbols
* `pytest test_kamodo.py::test_multi_arg_units`

### 2021-03-24 20:42:24.019445: clock-in

### 2021-03-24 18:52:07.659737: clock-out

* added test for multi argument unit composition

### 2021-03-24 18:51:37.128420: clock-in


```python
from sympy.abc import _clash
{'C': C,
 'O': O,
 'Q': Q,
 'N': N,
 'I': I,
 'E': E,
 'S': S,
 'beta': beta,
 'zeta': zeta,
 'gamma': gamma,
 'pi': pi}
```

### 2021-03-15 10:00:09.146598: clock-out


### 2021-03-15 09:59:36.883356: clock-in

### 2021-03-10 16:01:04.687207: clock-out

* adding citations
* pinning sympy for tests

### 2021-03-10 14:56:13.114546: clock-in: T-8m 

* updating sympy version

### 2021-03-10 12:21:29.440527: clock-out

* related projects

### 2021-03-10 12:20:49.749406: clock-in: T-80m 

### 2021-03-10 11:00:15.169767: clock-out

* fixed latex unit printing

### 2021-03-10 10:38:51.374290: clock-in

### 2021-03-09 20:00:17.427586: clock-out

* fixing latex unit rendering

### 2021-03-09 19:10:48.802378: clock-in

### 2021-03-09 12:35:42.893424: clock-out

* paper updates

### 2021-03-09 10:45:05.148443: clock-in

### 2021-03-08 20:17:26.525866: clock-out: T-70m 

### 2021-03-08 20:16:19.908822: clock-out


### 2021-03-08 17:43:02.224210: clock-in

### 2021-03-08 12:01:51.947367: clock-out


### 2021-03-08 11:23:30.640450: clock-in

### 2021-03-05 12:13:29.642736: clock-out

* started paper

### 2021-03-05 11:49:41.887072: clock-in

### 2021-02-16 20:33:17.491435: clock-out

* added zoom_test.html
* adding plot meta, made datetime test query deterministic

### 2021-02-16 20:03:14.752210: clock-in

### 2021-02-15 21:28:41.571728: clock-out


### 2021-02-15 20:41:15.499581: clock-in

### 2021-02-14 12:19:44.292293: clock-out: T-20m 


### 2021-02-14 10:13:19.121414: clock-in

### 2021-02-13 23:07:16.031488: clock-out


### 2021-02-13 22:48:04.844484: clock-in

### 2021-02-13 19:48:22.341859: clock-out

* js rangeslider widget

### 2021-02-13 17:49:49.306768: clock-in

### 2021-02-13 16:33:06.790805: clock-out


### 2021-02-13 13:24:50.034083: clock-in

### 2021-02-13 13:24:46.755317: clock-out: T-10m 


### 2021-02-13 13:01:49.691331: clock-in

### 2021-02-13 01:36:52.580893: clock-out

* developed slice widget for notebook

### 2021-02-12 22:29:50.872811: clock-in

### 2021-02-12 11:22:44.328681: clock-out


### 2021-02-12 10:31:25.170196: clock-in

* adding docs endpoint
### 2021-02-11 00:41:02.484844: clock-out

* fixed bug in kamodoAPI registering funtions twice

### 2021-02-10 23:11:21.144902: clock-in

### 2021-02-10 22:25:01.426986: clock-out

* function generator operations

### 2021-02-10 19:42:08.114657: clock-in

### 2021-02-09 13:11:09.281176: clock-out

* pushing apembroke/kamodo:0.1

### 2021-02-09 13:05:36.824475: clock-in

### 2021-02-09 13:05:10.883744: clock-out

* allowing function defaults with null args
* `flask.jsonify` should be used when returning from custom `get` methods

### 2021-02-09 12:50:49.247565: clock-in

### 2021-02-09 10:40:48.036338: clock-out

* need to include in `POST` method the ability to reference global models
* including user_model in kamodo.yaml, added delete method

### 2021-02-09 09:38:01.649942: clock-in

### 2021-02-09 01:21:34.140453: clock-out

* pinning sympy for api
* fixed bug in jsonifying user funcs

### 2021-02-08 23:05:27.925457: clock-in

### 2021-02-08 22:57:08.527566: clock-out


### 2021-02-08 22:48:42.433467: clock-in

### 2021-02-08 21:26:50.968524: clock-out

* need a way to post changes to variables
* added user user model endpoints
how do we determine default user model? - setting this in `kamodo.yaml`

got user endpoints working
* user models: `kamodo/api/usermodel`
* global models: `/api/modelA`


### 2021-02-08 17:24:04.082408: clock-in

### 2021-02-08 15:36:06.082620: clock-out

* registering custom fields
* need to deserialize as numpy arrays in combination with object_hook

## developer meeting
* emmpy - empirical modeling in python (tsygenenko, etc)
* pypluto
* ccmc docs have been updated

### 2021-02-08 13:07:56.346344: clock-in

### 2021-02-08 12:05:21.062862: clock-out

* user model endpoints

## user models
- `/kamodo/usermodelA/api`
- `/api/servermodelA`
- `/api/servermodelB`

### 2021-02-08 11:45:28.857313: clock-in

### 2021-02-08 10:38:31.202569: clock-out


### 2021-02-08 10:11:35.338486: clock-in

### 2021-02-06 11:37:06.027891: clock-out

* testing post method for user-defined expression

### 2021-02-06 11:01:27.615022: clock-in

### 2021-02-05 13:38:41.602116: clock-out

* prototyping user-defined kamodo objects
* added default forwarding for expressions

### 2021-02-05 11:10:12.407904: clock-in

* added data endpoint for cached function result
* installed requests for workflow
* kamodoAPI only registers units
* can add data to api to avoid initial call with defaults

### 2021-02-03 12:51:16.021583: clock-out

* serialized lambda generators

### 2021-02-03 12:11:53.187449: clock-in

### 2021-02-03 10:05:23.688588: clock-out

* serializing/deserializing generators - need to forge deserialized signatures? kamodofy?
* added flask.host and flask.port

### 2021-02-03 09:37:31.803474: clock-in

* test accessing multiple kamodo objects in same namespace
* fixed index serialization

### 2021-02-02 11:28:10.308687: clock-out

* improved serialization tests

### 2021-02-02 11:16:43.496588: clock-in

### 2021-02-02 09:37:36.981135: clock-out

* removing requests-mock
* more robust serialization/deserialization

### 2021-02-02 08:43:23.510270: clock-in

### 2021-02-01 19:27:49.355746: clock-out

* resolving serialization issues

### 2021-02-01 17:45:54.355466: clock-in

### 2021-02-01 08:59:03.221914: clock-out

* can use byte swapping and restore endianess https://numpy.org/doc/stable/user/basics.byteswapping.html
* jaweson https://github.com/someones/jaweson

### 2021-02-01 08:15:47.966367: clock-in

### 2021-01-30 14:13:25.401017: clock-out

trying different serialization methods

binary options:
* msgpack - binary serialization of numpy https://github.com/lebedov/msgpack-numpy https://github.com/msgpack/msgpack-javascript
* json and base64 encoding only https://stackoverflow.com/a/30698135
* bson `from bson.json_util import dumps, loads`

ascii options:
* pandas build_table_schema https://pandas.pydata.org/pandas-docs/version/0.21.0/generated/pandas.io.json.build_table_schema.html
* datapackage/table-schema https://github.com/frictionlessdata/datapackage-pipelines-pandas-driver#getting-started

### 2021-01-30 13:03:53.446801: clock-in

* fixed plot title latex
* added KamodoAPI class
* overriding models with config
* fixed datetime api example
* serialized time series, fixed default plot
* adding kTest to kamodo.yaml
* made test object have defaults

### 2021-01-21 23:35:17.186812: clock-out

* updated API.dockerfile
* added .dockerignore
* added function-specfic plot resource

### 2021-01-21 23:14:54.363159: clock-in

* added default str json
* added defaults output
* added kamodo-serve
* fixing error msg

### 2021-01-12 16:11:56.285579: clock-out

* checking for reserved names

### 2021-01-12 15:49:29.698018: clock-in

* cleaned up dimensioness latex

### 2020-12-21 13:52:51.879393: clock-out

* dimensionless composition passes

### 2020-12-21 13:32:17.655577: clock-in

* added unitless composition test
* fixed latex for kamodofied functions with no equations
* added arcseconds
* need to make sure kamodo objects remove signature when deleting items

### 2020-12-16 17:51:46.867654: clock-out

* added API.Dockerfile
* set ip for localhost in flask container, antlr
* merging util.py from psi/ensemble
* docker container

 | git merge master |	git rebase master
---- | ---- | ----
Keep changes from master	| --theirs	| --ours
Keep changes from feature |	--ours	| --theirs


### 2020-12-16 16:04:52.924287: clock-in

### 2020-12-16 14:56:58.695307: clock-out

* ccmc tagup
* `docker run -p 5000:5000 asherp/kamodo`
* adding flask dependency, running api on dockerfile startup
* adding dockerfile from nasa branch
* adding hapi

### 2020-12-16 13:10:31.733582: clock-in

* kamodo boolean operations?
* range queries?
* search?
* jupyter widgets - can prototype dashboards

### 2020-12-15 17:27:33.675261: clock-out

* merging test_utils

### 2020-12-15 16:32:25.927633: clock-in

### 2020-12-15 13:29:36.958990: clock-out

* adding coverage badge
* made repr_latex pass
* removed symbolic function
* test_repr_latex removing extra slash
* added decorator dependency
* graphviz - can be used to visualize the user's pipeline

```console
conda install -c conda-forge graphviz xorg-libxrender xorg-libxpm
pip install graphviz
```

### 2020-12-15 11:35:04.884012: clock-in

### 2020-12-14 18:29:33.877095: clock-out

* automated testing
* improved coverage

### 2020-12-14 17:41:49.324815: clock-in: T-35m 

### 2020-12-14 17:01:22.598481: clock-out


### 2020-12-14 16:44:16.293931: clock-in

### 2020-12-14 14:45:50.930880: clock-out

* fixed latex rendering

### 2020-12-14 13:35:34.718198: clock-in

### 2020-12-09 18:56:30.911547: clock-out

* inheriting from UserDict, latex tests
* article on alternatives to dict inheritance  https://treyhunner.com/2019/04/why-you-shouldnt-inherit-from-list-and-dict-in-python/


### 2020-12-09 16:46:32.320646: clock-in

### 2020-12-09 15:02:57.443215: clock-out

* improving coverage
* tests passing

### 2020-12-09 12:56:26.285722: clock-in

### 2020-12-08 16:29:15.942982: clock-out

* test_contains

### 2020-12-08 14:04:35.969497: clock-in

### 2020-12-07 18:07:09.300996: clock-out

* boosting code coverage

### 2020-12-07 15:47:35.229728: clock-in

### 2020-12-07 13:55:31.954661: clock-out

* added test for from_kamodo

### 2020-12-07 13:54:57.942787: clock-in

### 2020-12-07 12:02:40.382749: clock-out

* looking at unicode symbols

### 2020-12-07 11:04:24.226818: clock-in

### 2020-12-04 14:35:27.018180: clock-out

* increased coverage
Generating coverage report:

```bash
(kamodo)$ pytest --cov kamodo.kamodo --cov kamodo.util --cov plotting test_plotting.py test_kamodo.py util.py --cov-report html
```

### 2020-12-04 12:54:18.896202: clock-in

### 2020-12-04 12:48:07.464397: clock-out

* stress testing

### 2020-12-04 12:40:30.271441: clock-in

### 2020-12-03 16:33:17.672964: clock-out

* cleaned up conversion factors as fractions
* algebraic unit composition algebra

### 2020-12-03 15:45:57.247823: clock-in

### 2020-12-03 14:59:25.677055: clock-out


### 2020-12-03 14:31:39.023179: clock-in

### 2020-12-03 13:46:56.103902: clock-out

* adding algebraic unit composition tests
* removed resolve_unit in favor of get_expr_unit

### 2020-12-03 11:32:37.828682: clock-in

### 2020-12-02 17:11:43.336189: clock-out

* fixed bugs in base units
* passing tests

### 2020-12-02 15:09:11.720943: clock-in

### 2020-12-02 13:01:09.866258: clock-out

* isolated unit composition conversion issues

### 2020-12-02 10:52:59.764421: clock-in

### 2020-12-01 17:55:11.661674: clock-out

* resolving argument units

### 2020-12-01 16:29:52.875869: clock-in

### 2020-12-01 15:00:14.542857: clock-out

* bug fixing

### 2020-12-01 12:37:44.119942: clock-in

### 2020-11-30 19:05:54.863344: clock-out

* trying to resolve unit composition bugs

### 2020-11-30 17:00:58.216070: clock-in

### 2020-11-30 12:10:03.081932: clock-out

* finding bugs in unit composition

### 2020-11-30 10:56:13.021112: clock-in

### 2020-11-28 14:51:56.208988: clock-out

* clean up
* changing signature of unit composition

using a more intuitive signature

```python
kamodo['f(x[km],y[cm])[kg]'] = ...
```

### 2020-11-28 11:09:32.470142: clock-in

### 2020-11-27 15:41:15.004368: clock-out

* fixed unit conversion flip

### 2020-11-27 14:57:36.010904: clock-in

### 2020-11-27 14:21:05.568354: clock-out

* fixed composing with multiplied functions

### 2020-11-27 14:02:45.405415: clock-in

### 2020-11-26 14:57:57.717229: clock-out

* all kamodo tests passing

### 2020-11-26 14:08:27.224461: clock-in

### 2020-11-26 11:53:11.236857: clock-out

* raising NameError for bad unit conversion

### 2020-11-26 11:25:57.649820: clock-in

### 2020-11-25 18:57:00.223313: clock-out

* got unit composition to pass tests

### 2020-11-25 17:02:25.523707: clock-in

### 2020-11-25 14:28:39.430228: clock-out

* fixed vectorize composition

### 2020-11-25 13:36:07.998404: clock-in

### 2020-11-25 13:36:05.297272: clock-out: T-1h 


### 2020-11-25 12:19:16.144012: clock-in

### 2020-11-25 11:58:57.636064: clock-out


### 2020-11-25 11:42:44.748499: clock-in

### 2020-11-25 11:38:18.718909: clock-out


### 2020-11-25 10:55:45.291278: clock-in

### 2020-11-24 18:42:53.683721: clock-out

* unit composition

### 2020-11-24 16:51:00.759291: clock-in

### 2020-11-17 17:35:34.486919: clock-out


### 2020-11-17 16:06:38.651579: clock-in

### 2020-11-17 13:46:18.035319: clock-out

* passing functional unit test

### 2020-11-17 11:17:10.446484: clock-in

### 2020-11-16 18:34:07.584293: clock-out

* implementing unit functions

### 2020-11-16 16:59:30.697802: clock-in

### 2020-11-16 15:01:46.527485: clock-out


### 2020-11-16 13:11:21.254349: clock-in

### 2020-11-16 13:05:29.126719: clock-out: T-1h 


### 2020-11-16 11:12:27.592786: clock-in

### 2020-11-15 17:03:30.218928: clock-out


### 2020-11-15 15:17:08.285800: clock-in

### 2020-11-15 13:09:56.154786: clock-out


### 2020-11-15 11:52:34.355935: clock-in

### 2020-11-13 23:16:23.702830: clock-out


### 2020-11-13 22:37:29.395330: clock-in

### 2020-11-13 18:50:51.601791: clock-out


### 2020-11-13 17:01:21.677780: clock-in

### 2020-11-13 11:59:07.207912: clock-out

* overhauling units

### 2020-11-13 11:23:14.433529: clock-in

### 2020-11-12 15:25:13.259993: clock-out


### 2020-11-12 13:53:04.513345: clock-in

### 2020-11-12 11:40:29.694990: clock-out


### 2020-11-12 11:22:56.545479: clock-in

### 2020-11-11 17:42:08.055155: clock-out

* updated convert_to to raise errors

### 2020-11-11 16:46:46.677731: clock-in

### 2020-11-11 16:29:02.695711: clock-out


### 2020-11-11 16:10:38.700105: clock-in

### 2020-11-09 11:08:06.200728: clock-out

* fixing unit bugs

### 2020-11-09 10:00:15.420963: clock-in

* fixing unit conversion bugs

### 2020-11-03 11:17:20.736909: clock-out

* added dynamic function evaluation

### 2020-11-03 09:18:36.269712: clock-in

* radians
* why can't user access constants?

### 2020-10-27 23:26:44.585965: clock-out

* developing evaluate endpoint


`http://127.0.0.1:5000/api/mymodel/evaluate?variable=%27g=(f%2B1)**.5%27&x=[3,4,5]`

plus sign: `%2B`

### 2020-10-27 21:16:42.216480: clock-in

### 2020-10-27 15:30:20.449578: clock-out

* merging tests from Dhruv

### 2020-10-27 13:31:24.216328: clock-in

* code cleanup

### 2020-10-19 18:22:01.985135: clock-out

* all util.py unit tests pass
* spacing

* trying to fix collections warning

```bash
DeprecationWarning: Using or importing the ABCs from 'collections' instead of from 'collections.abc' is deprecated since Python 3.3,and in 3.9 it will stop working
```

### 2020-10-19 17:32:57.246898: clock-in

### 2020-10-14 23:22:57.047235: clock-out

* commenting util.py

### 2020-10-14 21:43:10.943127: clock-in

### 2020-10-14 12:31:50.544849: clock-out

* going through util.py

### 2020-10-14 11:29:55.853990: clock-in

### 2020-10-13 15:52:04.229990: clock-out

* cleaning up unit conversion code

### 2020-10-13 12:10:24.428950: clock-in

### 2020-09-02 16:03:51.164640: clock-out

* how to keep NASA readers and core from conflicting:
	- currently these are separate files so merges should be straight-forward
	- readers are all subclasses of Kamodo, so breaking changes should only be downstream
* modifications to core should be made with hourly to comply with NASA license
* could also rewrite kamodo core as a separate package with its own repo

### 2020-09-02 14:59:41.540200: clock-in

### 2020-08-12 12:54:08.590956: clock-out

* dev meeting

### 2020-08-12 11:59:10.698783: clock-in

### 2020-08-05 12:37:04.934436: clock-out

* dev meeting

### 2020-08-05 12:01:11.230431: clock-in

### 2020-07-29 12:58:25.948624: clock-out

* dev meeting
* modularity
	- make sure readers are independent
	- allow matplotlib without installing plotly

### 2020-07-29 12:02:08.702440: clock-in

### 2020-07-15 13:00:19.299408: clock-out

* dev meeting
* send blurb on cli to darren for gem

### 2020-07-15 12:31:56.828729: clock-in: T-15m 

* developer meeting
### 2020-07-08 12:31:42.478226: clock-out

* developer meeting
* attend GEM/AGU?

### 2020-07-08 11:59:12.669083: clock-in

### 2020-07-08 11:52:18.410732: clock-out


### 2020-07-08 11:50:08.658452: clock-in

### 2020-07-01 12:54:01.716646: clock-out

* developer meeting

### 2020-07-01 12:04:30.950750: clock-in

### 2020-06-24 12:54:16.183808: clock-out

* developer meeting

### 2020-06-24 12:16:33.594954: clock-in

### 2020-06-17 12:13:35.118654: clock-out


### 2020-06-17 12:12:50.788007: clock-in

### 2020-06-10 12:46:57.416636: clock-out

* developer meeting

### 2020-06-10 12:16:22.661881: clock-in

### 2020-06-03 13:01:51.004408: clock-out

* kamodo team meeting

### 2020-06-03 13:01:28.636263: clock-in: T-50m 

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

### 2020-05-27 14:38:38.650053: clock-out


### 2020-05-27 14:35:57.387894: clock-in

### 2020-05-27 12:57:54.286201: clock-out

* developer meeting

### 2020-05-27 11:38:16.287492: clock-in

### 2020-05-23 20:54:58.496791: clock-out

* fixed ordering bug in generator function evaluation

### 2020-05-23 20:01:46.123612: clock-in

* cleaning up skew contour carpet plots
* squeezing gridify output, added rvert lvert

### 2020-05-13 15:02:23.437971: clock-out

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

### 2020-05-13 09:06:37.493882: clock-in

### 2020-04-29 13:42:44.161370: clock-out

* flask server integration, api test working
* PYHC meeting

### 2020-04-29 12:14:17.372197: clock-in

### 2020-04-22 14:16:35.278555: clock-out: T-1h16m 

* More on flask integration from [plotly](https://dash.plotly.com/integrating-dash)


### 2020-04-22 11:56:00.335536: clock-in

### 2020-04-15 13:16:43.205737: clock-out: T-15m 


### 2020-04-15 11:57:50.421409: clock-in

### 2020-04-08 13:01:50.501409: clock-out

* meeting with developers

### 2020-04-08 12:18:16.752330: clock-in

### 2020-04-01 13:52:05.021374: clock-out

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

### 2020-04-01 12:08:40.151618: clock-in

### 2020-03-25 14:02:55.608031: clock-out

## kamodo meeting
* set deadlines for gui, april->may
* need a function that converts rho(x,y,z) -> rho(xvec)
* [dash-flask app tutorial](https://hackersandslackers.com/plotly-dash-with-flask/)
* tried running flask_restful with Dash, got `'Dash' object has no attribute 'handle_exception'`

### 2020-03-25 12:30:23.986927: clock-in

### 2020-03-25 12:20:20.625681: clock-out


### 2020-03-25 12:04:20.154862: clock-in

### 2020-03-18 21:37:24.225992: clock-out: T-5h 


### 2020-03-18 16:26:22.789889: clock-in

### 2020-03-18 13:03:13.421786: clock-out

* looking at flask rest api https://flask-restful.readthedocs.io/en/latest/quickstart.html#full-example

## developer meeting
* sscweb is positions also extrapolated into the future
* cdaweb is positions and data

### 2020-03-18 11:31:13.223644: clock-in

### 2020-03-11 20:25:31.313266: clock-out

* developer meeting

### 2020-03-11 20:25:18.435898: clock-in: T-45m 

### 2020-03-05 00:02:00.846556: clock-out

* prototyping mas

### 2020-03-04 22:01:13.683543: clock-in: T-30m 

### 2020-03-04 13:12:38.209902: clock-out

## developer meeting
* put website url in github description
* pyHC standards grading

https://github.com/heliophysicsPy/heliophysicsPy.github.io/blob/master/_pyhc_projects/pyhc_project_grading_guidelines.md

### 2020-03-04 11:52:51.855025: clock-in

### 2020-02-25 15:49:05.616947: clock-out

* addressing Liang's suggestions
* fixed deprecation warning from sympy>=1.5

### 2020-02-25 14:37:23.285614: clock-in

* check out partial functions
### 2020-02-19 13:03:32.424355: clock-out

* developer meeting

### 2020-02-19 12:22:27.823834: clock-in

### 2020-02-16 17:09:14.266657: clock-out: T-56h 

* submitted several feature issues

### 2020-02-14 09:56:13.530302: clock-in

### 2020-02-13 10:39:58.363603: clock-out

* iSWAT-COSPAR sessions
* end-to-end solutions
* community involvement
* visibility
* understanding interpolators

### 2020-02-13 09:09:49.084779: clock-in

### 2020-02-13 09:09:44.638954: clock-out: T-20h 


### 2020-02-12 09:37:08.349586: clock-in: T-30m 

### 2020-02-10 12:05:07.684918: clock-out

* fixed cli bug that prevented multiple plots from being saved

### 2020-02-10 12:03:49.344691: clock-in: T-1h 

### 2020-02-09 22:22:51.749709: clock-out

* gui and cli release
* fixed continuous reload bug

### 2020-02-09 20:13:13.525126: clock-in

### 2020-02-09 11:07:16.234918: clock-out


### 2020-02-09 10:30:57.198701: clock-in

### 2020-02-08 19:06:04.754325: clock-out

* gui improvements

### 2020-02-08 17:41:35.076690: clock-in

### 2020-02-08 15:54:29.107857: clock-out

* got reload config to work\!

### 2020-02-08 15:17:25.461178: clock-in

### 2020-02-08 14:26:00.479448: clock-out

* tried getting fully dynamic dash callbacks to work

### 2020-02-08 12:40:38.991214: clock-in

### 2020-02-07 19:58:47.600209: clock-out

* got interactiv configuration

### 2020-02-07 18:28:09.035029: clock-in

### 2020-02-05 18:28:53.256342: clock-out


### 2020-02-05 16:56:46.743551: clock-in

### 2020-02-05 16:50:08.675833: clock-out

* testing stateful storage update

### 2020-02-05 16:22:00.871936: clock-in

### 2020-02-05 14:14:54.356827: clock-out

* got clientside subplots to render
* meeting for iSWAT-COSPAR
* pushed code into NASA master
* need to update pypi version
* think about kamodo api
* docker container!

### 2020-02-05 11:23:08.377881: clock-in

### 2020-02-04 14:59:50.123294: clock-out

* got gui to load separate models and parameters

### 2020-02-04 13:22:32.543080: clock-in

### 2020-02-04 11:59:54.684191: clock-out

* reading up on dcc.store and clientside_callback

### 2020-02-04 10:39:40.567881: clock-in

### 2020-02-03 16:29:29.417205: clock-out

* added parameter checkboxes

### 2020-02-03 15:49:51.793947: clock-in

### 2020-02-03 12:36:08.959116: clock-out

* rendering equations through katex

### 2020-02-03 10:22:19.591903: clock-in

### 2020-01-31 19:35:30.915089: clock-out

* got range slider to work

### 2020-01-31 18:37:19.539544: clock-in

### 2020-01-31 17:56:33.727653: clock-out


### 2020-01-31 17:34:13.290486: clock-in

### 2020-01-31 17:19:16.956235: clock-out


### 2020-01-31 17:02:15.856951: clock-in

### 2020-01-31 16:52:07.773236: clock-out


### 2020-01-31 16:06:35.834953: clock-in

### 2020-01-31 15:11:15.638982: clock-out


### 2020-01-31 14:49:48.788349: clock-in

### 2020-01-31 13:27:55.782459: clock-out

* packaging pysatKamodo

### 2020-01-31 12:37:50.196753: clock-in

### 2020-01-31 11:19:55.778624: clock-out

* gui generating callbacks after layout is set

### 2020-01-31 10:35:50.687658: clock-in

### 2020-01-29 13:18:25.114169: clock-out

* meeting with Kamodo team for iSWAT-COSPAR prep

### 2020-01-29 11:59:08.216757: clock-in

* merging hourly.yaml
### 2020-01-29 11:19:29.007340: clock-out

* added link to NASA CCMC site
* cleaned up documentation site
* cleaning up docs

### 2020-01-29 09:50:29.055214: clock-in: T-20m 

### 2020-01-28 13:39:50.454679: clock-out

* answering support email

### 2020-01-28 13:30:23.612051: clock-in

### 2020-01-22 15:09:42.523706: clock-out


### 2020-01-22 15:09:36.423725: clock-in: T-40m 

* adding separate tabs for each model

### 2020-01-22 14:19:44.506406: clock-out

* working on gui layout

### 2020-01-22 12:13:09.639567: clock-in

### 2020-01-21 17:53:12.595685: clock-out

* added dynamic line plots to gui

### 2020-01-21 15:28:20.175086: clock-in

### 2020-01-21 13:37:38.575443: clock-out

* added composition to cli
* made compose use keyword arguments
* added compose function for multiple kamodo objects

### 2020-01-21 11:10:26.570097: clock-in

### 2020-01-20 15:23:39.237942: clock-out

* adding multiple models to cli
* setting up kamodo.yaml overrides
* Answered Lutz's issue

### 2020-01-20 11:21:44.788233: clock-in

### 2020-01-15 13:11:48.588900: clock-out

* testing hourly commit
* addressing issue from Lutz
* send agenda for meeting

### 2020-01-15 11:45:55.727665: clock-in

### 2019-12-20 13:50:20.467183: clock-out
* working on cli
* added ability to call kamodo from any work directory containing config.yaml
* added embedable plots: `model.plot_conf.output_type=div` saves the div to the file
* have markdown-include point to the div, so you can embed in your own docs!

### 2019-12-20 11:18:35.468801: clock-in

### 2019-12-19 12:17:51.613183: clock-out
* adding config file description to docs
* Looking at `Config search path`. 
* hydra plugin configs are found automatically once they are installed
* how do we use the cwd to configure? https://github.com/facebookresearch/hydra/issues/274
* docusaurus looks interesting

### 2019-12-19 11:26:13.027723: clock-in

### 2019-12-09 13:53:18.142633: clock-out
* forgot to include cli notebook

### 2019-12-09 13:52:30.292970: clock-in

### 2019-12-09 13:50:43.574951: clock-out
* working on cli

### 2019-12-09 10:05:45.828781: clock-in

### 2019-12-05 21:55:15.506661: clock-out
* looking at Y-combinator for recursive anonymous functions

### 2019-12-05 21:54:44.630995: clock-in: T-90m

### 2019-11-27 12:42:29.395096: clock-out
* fixing cdf arrays

### 2019-11-27 12:06:51.013991: clock-in

### 2019-11-22 09:57:46.571462: clock-out

* regex matching https://stackoverflow.com/questions/1687620/regex-match-everything-but-specific-pattern

### 2019-11-22 09:09:13.743347: clock-in

### 2019-11-21 12:09:12.103191: clock-out

### 2019-11-21 11:01:35.147547: clock-in

### 2019-11-21 10:27:47.314091: clock-out
* cdflib: switching to multiindex for all dependencies

### 2019-11-21 09:50:46.684607: clock-in

### 2019-11-20 16:27:17.501485: clock-out
* moved docs/notebooks/kameleon/kameleon_gateway.py into kamodo/readers
* pushed recent work
* created pull request

### 2019-11-20 15:30:59.817396: clock-in

### 2019-11-20 12:42:32.541944: clock-out

* on forwarding defaults: consider leveraging a function's .data during composition.
this would allow each downstream function to have automatic defaults!
* need to push these changes

### 2019-11-20 12:06:16.296442: clock-in

### 2019-11-20 11:58:22.508485: clock-out
* Kamodofied cdflib!

### 2019-11-20 09:35:49.895291: clock-in

### 2019-11-13 09:40:03.974566: clock-out

### 2019-11-13 08:58:52.993108: clock-in

### 2019-11-12 20:37:52.561965: clock-out

### 2019-11-12 19:17:10.399952: clock-in

### 2019-11-12 18:55:01.500196: clock-out
* command line

### 2019-11-12 18:37:46.018928: clock-in

### 2019-11-11 10:20:48.198321: clock-out
* set up command-line plotting

### 2019-11-11 09:35:53.331821: clock-in

### 2019-11-08 16:51:57.754553: clock-out
* begin work on command line interface
* trying facebook's hydra cli architecture
* `pip install hydra-core --upgrade`

### 2019-11-08 15:53:02.859151: clock-in

### 2019-10-18 16:47:23.366813: clock-out
* looking at inverse mapping
* Kamodofied pytiegcm
* need to finish gridifying inverse mapping
* consider adding indexing option to gridify

### 2019-10-18 11:58:47.333047: clock-in

### 2019-10-17 17:54:11.880758: clock-out
* kamodofying pyTIEGCM

### 2019-10-17 17:34:05.273266: clock-in

### 2019-10-17 15:54:57.284128: clock-out
* Kamodofied pyTIEGCM

### 2019-10-17 13:20:49.585114: clock-in

### 2019-09-30 12:34:55.407880: clock-out: added kamodofied kameleon example
* Kameleon-kamodo bridge development
* Added kamodofied kameleon object

### 2019-09-30 10:00:33.528255: clock-in

### 2019-09-11 13:32:56.778768: clock-out
* Created Fieldline Tutorial
* made griddify return in `xy` indexing to work with map-to-plane
* added scatter plot to available plot_types
* Finished kamodofying_models tutorial
* Need to check if plotting map-to-plane works for each axis

### 2019-09-11 08:19:33.545586: clock-in

### 2019-09-10 14:49:49.344817: clock-out
* fixed cone colors in plotting

### 2019-09-10 11:39:05.825143: clock-in

### 2019-09-10 10:43:05.346628: clock-out
* removing nans from solver output in favor of MultiIndex seeds

### 2019-09-10 09:55:34.659931: clock-in

### 2019-09-09 13:43:18.711955: clock-out
* Added pandas i/o for 3d line and vector plots

### 2019-09-09 11:07:35.809509: clock-in

### 2019-09-09 10:52:08.229479: clock-out
* Added pointlike decorator

### 2019-09-09 09:37:01.946111: clock-in

### 2019-09-06 15:26:26.921477: clock-out
* Looking into LFM wrapper
* https://wiki.ucar.edu/display/LTR/pyLTR seems to have been written for python 2

### 2019-09-06 15:20:58.638859: clock-in

### 2019-09-05 17:04:36.774012: clock-out
* solver decorator
* dipole field test
* stopping integration at boundary
* can choose resolution of knots

### 2019-09-05 11:09:34.195206: clock-in

### 2019-09-04 15:29:02.593470: clock-out
* got fieldline tracer working (ivp solver)

### 2019-09-04 14:47:48.601879: clock-in

### 2019-09-04 10:57:19.685814: clock-out

### 2019-09-04 10:53:23.666033: clock-in

### 2019-09-04 10:00:42.819821: clock-out
* looking at complex parameters for streamlines

### 2019-09-04 09:42:08.493375: clock-in

### 2019-09-03 18:33:35.737108: clock-out
* working on streamlines

### 2019-09-03 16:44:51.693606: clock-in

### 2019-09-03 14:09:22.239226: clock-out: T-6m
* Deduct 6m
* Handled pandas output for plotting 3d vector fields
* Fixing integration symbol for solver
* Fixed griddable vector fields

### 2019-09-03 13:07:54.277159: clock-in

### 2019-08-30 10:49:09.145324: clock-out

* Issue when calling util.parse_latex:

> ImportError: LaTeX parsing requires the antlr4 python package, provided by pip (antlr4-python2-runtime or antlr4-python3-runtime) or conda (antlr-python-runtime)

* reconsider `antlr4` dependency

### 2019-08-30 10:47:31.579538: clock-in

### 2019-08-28 16:30:55.501447: clock-out: added gridify decorator
* Expanding kamodofication tutorial to support more variables
* Looking at Numpy's generalized universal function API - would formalize our mapping to plotting routines
* Add example of point interpolation for trajectory

### 2019-08-28 10:24:56.645799: clock-in

### 2019-08-16 16:27:22.801411: clock-out: cleaned up imports, added Kamodofication example, bug fixes
* testing install without psutil
* creating conda package
* cleaned up imports
* added click dependency
* updated pip

### 2019-08-16 14:25:59.427481: clock-in

### 2019-08-14 16:18:12.261596: clock-out
* finishing kamodofied models example

* Had issues installing locally with `pip install .` Need to use `pip install -e .` instead.

### 2019-08-14 10:37:19.291831: clock-in

### 2019-08-12 12:14:17.649049: clock-out
* Adding example notebook on kamodofication

### 2019-08-12 10:27:35.850505: clock-in

### 2019-07-24 16:44:02.434437: clock-out
* Moving code into new repo
* Creating new kamodo-user environment from clean miniconda

### 2019-07-24 14:33:41.015853: clock-in

