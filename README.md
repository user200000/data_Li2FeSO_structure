# "Anion-polarisation--directed short-range-order" Analytical Workflow

Authors:
- Samuel W. Coles (Wrote workflow)
- Viktoria Falkowski
- Harry S. Geddes
- Gabriel E. Pérez
- Samuel G. Booth
- Alexander G. Squires 
- Conn O'Rourke
- Kit McColl
- Andrew L. Goodwin
- Serena A. Cussen
- Simon J. Clarke
- M. Saiful Islam
- Benjamin J. Morgan

## Summary


This repository contains the analytical workflow for our preprint entitled "Anion-polarisation--directed short-range-order in antiperovskite Li2FeSO” [LINK GOES HERE]. The workflow consists of three main components. Fitting, Enumeration and Monte Carlo (MC), which are each applied to coulombic and cluster expansion system. A random arrangement of cations is also generated for comparison with the monte carlo. There are three main subsections of the workflow fitting, mc, and enumeration each will be described briefly there is a workflow diagram at the base of the page. A conda yaml file is provided describing the exact environments in which calculations were initially run.

## Fitting

Fitting and all subsequent calculations are performed using the [icet library](https://icet.materialsmodeling.org). Data from DFT calculations are input as jsons and fit using the python scripts. The output from this step is cluster expansion object for DFT energies. In parallel a cluster expansion is fit using Ewald energies in the coulomb tab.

## MC

Monte carlo simulations are performed for 8x8x8 supercells using the two fitted cluster expansions using the mchammer library installed as part of icet. The monte carlo simulation outputs mc data container objects which are then processed by polyhedral.py using the [polyhedral analysis library] (https://polyhedral-analysis.readthedocs.io/en/latest/) to octahedral population figures.

## Enumeration

The final stage, enumeration, focusses on the direct calculation of the partition function. In enumeration.py the different supercells of a 2x2x2 expansion of the Li2FeSO primitive are obtained using [bsym](https://joss.theoj.org/papers/10.21105/joss.00370), the degeneracy of structures is obtained using this process. The energies of these structures are calculated along with key structural parameters used to obtain order parameters. In two subsequent scripts violins.py and paramaters_and_cdos.py we obtain Fig. 3 and Fig. 7 directly. paramaters_and_cdos.py handles the calculation of temperature dependent order parameters and the thermal population density of states at 1025K for both cluster expansions.

## Workflow scheme



![](./workflow.png)


An associated dataset will be made available with DFT calculations on the bath university archive.
