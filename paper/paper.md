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

Kamodo is a functional api for scientific models and data. In Kamodo,
all scientific resources are registered as symbolic fields which are typically mapped to
underlying model and data interpolators. This allows many common problems, such
as field line integration, coordinate transformation, and recontextualization,
to be posed as function compositions. Kamodo employs a transparent unit conversion scheme
that mimics hand-written expressions: units are declared on the left hand side of user
expressions via bracket notation and conversion factors are automatically inserted on
the right hand side. Kamodo includes quick-look graphics, dashboards, and a LaTeX interface,
and is amenable to containerization and cloud hosting. While Kamodo was designed
to solve the interperational challenges of the space weather industry, it is general
enough to be applied in other contexts.

# Statement of need

Space weather models and data employ a wide variety of specialized data formats,
data structures, and interfaces. Often these are suitable for
the needs of individual research groups, which have sophisticated pipelines tailored
for their field of study. However, this specialization poses an impediment
for research to operations, space weather forecasting, research, and education.
Kamodo is a symbolic abstraction layer that utilizes model and data interpolators
to fascilitates downstream visualization, recontextualization and
application development. Kamodo is built on Python libaries Sympy, Plotly, but
the general approach may be adapted to other langvisualization paradigms
Many common problems in science discovery. may be posed as function compositions.


# Unit System

Kamodo's unit system differs dramatically from many other powerful packages used
in space weather such as Astropy, for two reasons: First, Sympy includes its own
unit system, so reliance on a parallel unit system would be redundant.
Second, Kamodo's use of expressions allows units to be kept separately from the 
data types used by functions. The only requirement is that the types returned
support algebraic manipulation (e.g. support multiplication, addition, etc). 

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
