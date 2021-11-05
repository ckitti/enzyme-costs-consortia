# Use of SteadyCom for synthetic microbial consortia analysis

  * Repository Information

This repository contains the script used for the evaluation of a hypothetical four-member consortia based on auxotrophic strains of _E. coli_ using an updated GEM(*i*ML1515) and figures that accompany the technical report for this method for the course Advanced Technologies for Bioscience.

Additionally, the use of SteadyCom is also adapted to accommodate the use of the enzyme-constrained models ([ecModels](https://github.com/SysBioChalmers/ecmodels)) to investigate the influence of protein cost on community behavior and composition. A pipeline script is also available to study pairwise setups based on auxotrophic strains investigated by [Mee et al., 2014](https://doi.org/10.1073/pnas.1405641111).

This repository is administered by Cheewin Kittikunapong, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Dependencies

### Recommended Software:
* A functional Matlab installation (MATLAB 7.3 or higher)
* [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/)
* libSBML MATLAB API ([version 5.16.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended)
* [IBM ILOG Cplex ](https://opencobra.github.io/cobratoolbox/latest/installation.html#ibm-ilog-cplex)
* While not necessary, some supplemental analysis is done using the [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)
