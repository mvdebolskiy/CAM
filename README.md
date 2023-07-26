# CAM-Nor: The NorESM version of the Community Atmosphere Model (CAM)

CAM code is stored in this repository on branches other than main.  The details are explained below.

Please see the [NorESM documentation](https://noresm-docs.readthedocs.io/en/noresm2/) for complete documentation on NorESM, getting started with git and how to contribute to CAM-Nor's development.

This repository has the following CAM-Nor branches:

* **cam_cesm2_1_rel_05-Nor** - contains the CAM6-Nor code, based on CESM2.1 CAM code. Contains aerotab and further changes as described by [Seland et al. 2020](https://gmd.copernicus.org/articles/13/6165/2020/), with tags

  - cam_cesm2_1_rel_05-Nor_v1.0.0 -- CAM6-Nor version (NorESM tag [release-noresm2.0.1](https://noresm-docs.readthedocs.io/en/noresm2/access/releases_noresm20.html#noresm2-0-1)), initial preliminary release of CMIP6 version, 18 March 2020

  - cam_cesm2_1_rel_05-Nor_v1.0.1 -- CAM6-Nor version, not used in a NorESM tag, 11 May 2020

  - cam_cesm2_1_rel_05-Nor_v1.0.2 -- CAM6-Nor version (NorESM tag [release-noresm2.0.2](https://noresm-docs.readthedocs.io/en/noresm2/access/releases_noresm20.html#noresm2-0-2)), CMIP6 bit-identical version (if used on fram, and as MM), 17 July 2020

  - cam_cesm2_1_rel_05-Nor_v1.0.3 -- CAM6-Nor version (NorESM tag [release-noresm2.0.3](https://noresm-docs.readthedocs.io/en/noresm2/access/releases_noresm20.html#noresm2-0-3)), ??, 1 April 2021

  - cam_cesm2_1_rel_05-Nor_v1.0.4 -- CAM6-Nor version (NorESM tag [release-noresm2.0.4](https://noresm-docs.readthedocs.io/en/noresm2/access/releases_noresm20.html#noresm2-0-4)), ??, 8 May 2021

  - cam_cesm2_1_rel_05-Nor_v1.0.5 -- CAM6-Nor version (NorESM tags [release-noresm2.0.5](https://noresm-docs.readthedocs.io/en/noresm2/access/releases_noresm20.html#noresm2-0-5) and [release-noresm2.0.6](https://noresm-docs.readthedocs.io/en/noresm2/access/releases_noresm20.html#noresm2-0-6)), ??, 9 December 2022

* **cam_cesm2_1_rel_05-Nor-SpAer** - based on cam_cesm2_1_rel_05-Nor, about to contain simple plume aerosol model code

* **cam_cesm2_1_rel_05-Nor_v1.0.2_keyClim-withoutanthraer** - based on cam_cesm2_1_rel_05-Nor (branched off at tag cam_cesm2_1_rel_05-Nor_v1.0.2) : related to version used in project KeyClim for simulation without anthropogenic aerosol (only natural aerosol)

* **cam-Nor-dev** - based on cam_cesm2_1_rel_05-Nor (branched off at tag cam_cesm2_1_rel_05-Nor_v1.0.2) : post-CMIP6 development version of CAM6-Nor (including bugfixes, new developments, answer-changes modifications, ...)

* **feature-hamocc-vsls** - based on cam_cesm2_1_rel_05-Nor (branched off at tag cam_cesm2_1_rel_05-Nor_v1.0.2) : related to version developed in project KeyClim for simulation of Very Short Lived Spiecies (VSLS) in ocean and atmosphere.

* **noresm_develop** - main development branch for future release versions of NorESM

* **main** - contains this README, basic guidelines, and some special GitHub scripts.

### CAM-Nor and NorESM documentation https://noresm-docs.readthedocs.io/en/noresm2/

## How to checkout and use CAM:

The instructions below assume you have cloned this repository and are in the repository directory. For example:
```
git clone https://github.com/NorESMhub/CAM
cd CAM
```

### To run CAM compatible with the NorESM2 release:
```
git checkout cam_cesm2_1_rel_05-Nor_v1.0.5
./manage_externals/checkout_externals
```
### To view the release branch in github, go to the "Branch:main" pulldown menu and select cam_cesm2_1_rel_05-Nor

### To use unsupported CAM **development** code:

## NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](https://www.cgd.ucar.edu/sections/cseg/policies).
```
git checkout noresm_v4_cam6_3_112
./manage_externals/checkout_externals
```
### CAM6-Nor Documentation - https://noresm-docs.readthedocs.io/en/noresm2/model-description/atm_model.html

### CAM6 namelist settings - https://docs.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html
