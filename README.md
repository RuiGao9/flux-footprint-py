[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19076529.svg)](https://doi.org/10.5281/zenodo.19076529)
![Stars](https://img.shields.io/github/stars/RuiGao9/flux-footprint-py?style=social)<br>

## flux-footprint-py: A Toolkit for Kljun-based Footprint Modeling
Calculates spatial footprints for water-energy-carbon flux research using the Kljun et al. (2015) model, supporting monitoring of tower fetch at different temporal scales.A Python workflow for calculating and aggregating spatial footprints of water, energy, and carbon fluxes. Implements the Kljun (2015) model to monitor tower fetch across various temporal scales.

## FFP input
This section can refer to the document [A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP)](https://footprint.kljun.net/downloads/v1.42/FFP_readme_Python.pdf)<br>
All inputs as scalars  
- $z_m$ = Measurement height above displacement height (i.e. z-d) [m] 
- $z_0$ = Roughness length [m] - enter [None] if not known 
- $u_{mean}$ = Mean wind speed at zm [ms-1] - enter [None] if not known 
- $h$ = Boundary layer height [m] 
- $ol$ = Obukhov length [m] 
- $\sigma_v$ = Standard deviation of lateral velocity fluctuations [ms-1] 
- $u^*$ = Friction velocity [ms-1]  
Note: Either $z_0$ or $u_{mean}$ is required. If both are given, $z_0$ is selected to calculate the foorprint.<br>

Optional input
- $wind_{dir}$ = Wind direction in degrees (of 360) for rotation of the footprint 
- $rs$ = Percentage of source area, i.e. a value between 10% and 90%.  Can be either a single value (e.g., "80") or an array of increasing percentage values  (e.g., "[10:10:80]"). Expressed either in percentages ("80") or in fractions of 1 ("0.8"). Default is [10:10:80]. Set to "NaN" for no output of percentages. 
- $nx$ = Integer scalar defining the number of grid elements of the scaled footprint. Large nx results in higher spatial resolution and higher computing time. Default is 1000, nx must be >=600. 
- $rslayer$ = Calculate footprint even if zm within roughness sublayer: set rslayer = 1. Note that this only gives a rough estimate of the foorprint as the model is not valid within the roughness sublayer. Default is 0 (i.e. no figure).
- $crop$ = Crop output area to size of the 80% footprint or the largest r given if crop = 1
- $fig$ = Plot an example figure of the resulting footprint on the screen: set fig = 1. Default is 0 (i.e. no figure).

## Installation
```
pip install "git+https://github.com/RuiGao9/flux-footprint-py.git" 
```

## How to run this model
`run_fpt_model.ipynb` is the main program you can edit
- Call input (e.g., `example_meteo_days.txt` under `example_data` folder)
- You can define the folder path in the first section
- You can edit parameters in the second section
- Call the footprint model under `flux_footprint` folder in the third section and generate the footprint raster data

## Reference
<a name="ref1"></a>
1. Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P.: A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP), Geosci. Model Dev., 8, 3695–3713, https://doi.org/10.5194/gmd-8-3695-2015, 2015.<br>
2. Nassar A, Torres-Rua A, Kustas W, Nieto H, McKee M, Hipps L, Alfieri J, Prueger J, Alsina MM, McKee L, Coopmans C, Sanchez L, Dokoozlian N. To What Extend Does the Eddy Covariance Footprint Cutoff Influence the Estimation of Surface Energy Fluxes Using Two Source Energy Balance Model and High-Resolution Imagery in Commercial Vineyards? Proc SPIE Int Soc Opt Eng. 2020 Apr-May;11414:114140G. https://doi.org/10.1117/12.2558777. Epub 2020 May 26. PMID: 33758459; PMCID: PMC7982303.<br>
3. Gao, R., A. Nassar, A. F. Torres-Rua, L. Hipps, M. Aboutalebi, W. A. White, M. Anderson, W. P. Kustas, M. M. Alsina, J. Alfieri, N. Dokoozlian, F. Gao, H. Nieto, L. McKee, J. H. Prueger, L. Sanchez, A. J. Mcelrone, N. B. Ortiz, I. Gowing, C. Coopmans (2021). Footprint area generating based on eddy covariance records, HydroShare, https://doi.org/10.4211/hs.9118e2c1034e40e4ba4721cd17702f70

## How to cite this work
Gao, R., Khan, M., & Viers, J. (2026). flux-footprint-py: A Python Toolkit for the Implementation of the Kljun Footprint Model (initials). Zenodo. https://doi.org/10.5281/zenodo.19076529

## Repository update information
- Creation date: 2026-03-17
- Last update: 2026-03-17

## Contact inforamtion if issues were found
Rui Gao: <br>Rui.Ray.Gao@gmail.com or RuiGao@ucmerced.edu
