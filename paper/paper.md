---
title: 'Kamodo: A functional api for space weather models and data'
tags:
  - Python
  - plasma physics
  - space weather
authors:
  - name: Asher Pembroke
    affiliation: 1
  - name: Darren DeZeeuw
    affiliation: 2, 3
  - name: Lutz Rastaetter
    affiliation: 2
  - name: Rebecca Ringuette
    affiliation: 2
affiliations:
 - name: Predictive Science, Inc.
   index: 1
 - name: Community Coordinated Modeling Center, NASA GSFC
   index: 2
 - name: University of Michigan
   index: 3
date: Mar 5, 2021
bibliography: paper.bib
---

# Summary

Kamodo is a functional programing interface for scientific models and data. In Kamodo, all scientific resources are registered as symbolic fields which are mapped to model and data interpolators or algebraic expressions. Kamodo performs function composition and employs a unit conversion system that mimics hand-written notation: units are declared in bracket notation and conversion factors are automatically inserted into user expressions. Kamodo includes a LaTeX interface, automated plots, and a browser-based dashboard interface suitable for interactive data exploration. Kamodo's json API provides context-dependent queries and allows compositions of models and data hosted in separate docker containers. Kamodo is built primarily on sympy [@10.7717/peerj-cs.103] and plotly [@plotly]. While Kamodo was designed to solve the cross-displinary challenges of the space weather community, it is general enough to be applied in other fields of study.

# Statement of need

Space weather models and data employ a wide variety of specialized formats, data structures, and interfaces tailored for the needs of domain experts. However, this specialization is also an impediment to cross-displinary research. For example, data-model comparisons often require knowledge of multiple data structures and observational data formats. Even when mature APIs are available, proficiency in programing languages such as python is necessary before progress may be made. This further complicates the transition from research to operations in space weather forecasting and mitigation, where many disparate data sources and models must be presented together in a clear and actionable manner. Such complexity represents a high barrier to entry when introducing the field of space weather to newcomers at space weather workshops, where much of the student's time is spent installing prerequisite software. Several attempts have been made to unify all existing space weather resources around common standards, but have met with limited success. 

Kamodo all but eliminates the barrier to entry for space weather resources by exposing all scientifically relavent parameters in a functional manner. Many problems in space weather analysis, such as field line tracing, coordinate transformation, and interpolation, may be posed in terms of function composition, kamodo is and ideal tool in the scientist's workflow. Kamodo builds on existing standards and APIs and does not require programing expertise on the part of end user. Kamodo is expressive enough to meet the needs of most scientists, educators, and space weather forecasters, and Kamodo containers enable a rapidly growing ecosystem of interoperable space weather resources. 


# Related Projects

Kamodo is designed for compability with python-in-heliosphysics [@ware_alexandria_2019_2537188] packages, such as PlasmaPy [@plasmapy_community_2020_4313063] and PySat [@Stoneback2018], [@pysat200]. This is accomplished through Kamodo subclasses, which are responsible for registering each scientifically relevant variable with an interpolating function. Metadata describing the function's units and other supporting documentation (citation, latex formating, etc), may be provisioned by way of the `@kamodofy` decorator.

Kamodo's unit system is built on SymPy [@10.7717/peerj-cs.103] and shares many of the unit conversion capabilities of `Astropy` [@astropy] with two key differences: first, Kamodo uses an explicit unit conversion system, where units are declared during function registration and appropriate conversion factors are automatically inserted on the right-hand-side of final expressions, which permits back-of-the-envelope validation. Second, units are treated as function metadata, so the types returned by functions need only support algebraic manipulation (Numpy, Pandas, etc). Output from kamodo-registered functions may still be cast into other unit systems that require a type, such as Astropy [@astropy], Pint [@pint], etc.

Kamodo can mimic some of the capabilities of raw data APIs such as HAPI. As with other PyHC projects, the goal is not to replace existing APIs, but rather to extend their capabilities. For example, Kamodo's API support purely functional data access, where `GET` requests can specify positions or times for which interpolated values should be returned. In addition, Kamodo `POST` requests may be used to register new functions on the server which are compositions of previously defined variables and with custom units.

Kamodo container services may be built on other containerized offerings. Containerization allows dependency conflicts to be avoided through isolated install environments. Kamodo extends the capabilities of space weather resource containers by allowing them to be composed together via the KamodoClient, which acts as a proxy for the containerized resource running the KamodoAPI. 


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

# Figures


# Acknowledgements



# References