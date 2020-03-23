# Function Space Optimization: A symbolic regression method for estimating parameter transfer functions for hydrological models

<p align="center">
  <img width="460" src="https://github.com/MoritzFeigl/FSO_paper/blob/master/FSO.png">
</p>

Accompanying code for the publication "Function Space Optimization: A symbolic regression method for estimating parameter transfer functions for hydrological models"

The code in this repository was used to produce all results and figures in our manuscript. The data for geo-physical properties is available under https://doi.org/10.5281/zenodo.3676053 and the discharge data used in this study is available at https://www.ehyd.gv.at. The meteorological data from the INCA dataset cannot be made public, because the rights belong to the Zentralanstalt für Meteorologie und Geodynamik (ZAMG). It can be acquired from https://www.zamg.ac.at.

### Content of the repository

- `Functions` All functions used in the Paper code. This mainly includes functions for the context free grammar implementation, FSO VAE training and generator functions, d-GR4J optimization routines and additional helper functions. 
- `Paper code` The scripts used to produce the results of the publication.


### Additional informations

All scripts in this repository use tensorflow 1. We are now working on a implementation of FSO as an R & Python package which are using tensorflow 2 and will therefore not update the here presented code base.

## License

[Apache License 2.0](https://github.com/MoritzFeigl/FSO_paper/blob/master/LICENSE)
