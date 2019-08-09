# About Kamodo

Kamodo has been under development at the Community Coordinated Modeling Center (CCMC), NASA GSFC since May, 2018. The CCMC supports the space weather community by providing software and services guided by domain experts in a variety of heliophysic science domains.

Kamodo supports the goals of the CCMC by:

* Bringing together models and data into a single high-level mathematical framework
* Allows scientists and educators to work with complex space weather models and data with little or no coding experience
* Provides an easy-to-extend framework for developers.

### Kameleon legacy

Kamodo shares some similarities with its predecessor, the CCMC's Kameleon Software Suite, insofaras it provides a unified API for space weather models. However, Kamodo gets there through a very different means: by leveraging cutting-edge python projects from both the heliophysics community (sunpy, spacepy, etc.) as well as more general mathematical frameworks like sympy. This allows Kamodo to be much more broad in its application, able to handle arbirtary scientific data and physics-based models. At the same time, by  building on the tools provided by model and data providers, Kamodo inherits the high performance necessary for data analysis. We felt that due to the large departure in both design and scope from Kameleon, it was necessary to launch Kamodo as a separate project under a different moniker.

## Design philosophy

Primary Design considerations

* Open Source (Apache 2.0)
* Should be format-, model-, data-agnostic
* Should support all types of users (non-coders, devs, modelers). *Anyone can cook!*
* Carrot approach to metadata (useful, but not mandatory)
