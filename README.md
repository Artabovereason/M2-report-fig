
# Code used for the internship

The present scripts were developped during my M2 internship at Institut de Physique et de Chimie des Matériaux de Strasbourg (IPCMS) in the Quantum dynamics of nano-objects team (QDYNO). The aim of those codes is to compute the relevent quantity for either the classical capture model or the semi-classical one.

## Usage
* The script `report_photoionization.py` will auto-compile the Fortran code available on the repository `Code/` and read the output to then make a plot of the photoionization at a given energy given the first and second quantum number. This script is only used when we want to see the photoionization. In the other code where we need the photoionization, this Fortran program is called and then interpolate. The photoionization code was developped my K. Levêque.

* The script `report_Langevin.py` is simply the numerical plot of the Langevin cross-section.  

* The script `report_dipolarcoupling.py` read data on a numerised table (`dipole_strengths_data.xlsx`) that is from D. E. Ramaker and J. M. Peek, 1973.

* The script `true.py` generate the cross-section value for multiple principal quantum number. Be aware that this code can be lengthy calculation wise, that's why on this I only made 31 points per entry-state, but one could optimize it to gain precision. To plot it alongside the experimental data (`H(1s)+H(2s).csv`, `H(1s)+H(3s).csv` and `H(1s)+H(4s).csv`) available, `plot_quantities.py` is used.
## Authors

- Pierre Guichard

## License
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
