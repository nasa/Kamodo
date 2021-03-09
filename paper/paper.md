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


# Summary

Kamodo is a functional programing interface for scientific models and data. In Kamodo, all scientific resources are registered as symbolic fields which are mapped to algebraic expressions or to model and data interpolators. Kamodo's functional design allows many common problems, such as field line integration and coordinate transformation, to be posed in terms of function compositions familiar to scientists. Kamodo employs a unit conversion system that mimics hand-written expressions: units are declared via bracket notation and conversion factors are automatically inserted on the right hand side of user expressions. Kamodo includes a LaTeX interface, automated plots via plotly, and a browser-based dashboard interface suitable for interactive data exploration. Kamodo's json API provides context-dependent queries and allows compositions of models and data hosted in separate containers. While Kamodo was designed to solve the cross-displinary challenges of the space weather community, it is general enough to be applied in other fields of study.

# Statement of need

Space weather models and data employ a wide variety of specialized data formats, data structures, and interfaces tailored for many different domains of the heliosphere and specialized to meet the needs of individual research groups. However, this specialization is also an impediment to cross-displinary research. For example, data-model comparisons require knowledge of both model and observational data formats. Even when a mature API is available, proficiency in programing languages such as python is required before progress can be made. A similar difficulty arises in the transition from research to operations (R2O) in space weather forecasting, where many disparate data sources and models must be presented together in a clear and actionable manner. Such low-level complexity represents a high barrier to entry when introducing the field of space weather to newcomers during space weather workshops, where much of the student's time is spent installing prequisite software. Several attempts have been made to unify all existing space weather resources around a common data/metadata standard, but have met with limited success. Similarly, many object-oriented APIs have attempted to support novel data structures outside of the scope of their design. We propose a solution that builds on existing standards and APIs without imposing a new standard and without requiring programing expertise on the part of end users, yet is expressive enough to meet the needs of many scientists, educators, and space weather forecasters.


# CCMC

The Community Coordinated Modeling Center (CCMC) at NASA, GSFC provides computational resources, research-focused services, and domain expertise for a large number of space weather models.  

# Unit System

Kamodo's unit system differs dramatically from many other powerful packages used in space weather such as Astropy, for two reasons: First, Sympy includes its own unit system, so reliance on a parallel unit system would be redundant. Second, Kamodo's use of expressions allows units to be kept separately from the  data types used by functions. The only requirement is that the types returned support algebraic manipulation (e.g. support multiplication, addition, etc). 

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

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
