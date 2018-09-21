# Data from article Stocker et al. (2018b) *Nature Geosci.*

The datasets provided here include:

- Site-level FLUXNET GPP estimates and GPP model results from the P-model (Wang et al., 2017) and MODIS MOD17A2H (Running et al., 2004)
- Model outputs from global simulations with the P-model (Wang et al., 2017) as implemented for the study by Stocker et al. (2018b)

This data may be used to partly reproduce results presented in Stocker et al. (2018b) *Nature Geosci.*. "Partly" because we used data for our analysis that was not open access but was confidentially shared with us. This includes remote sensing-based GPP estimates from the BESS and VPM models. MODIS GPP was downloaded from an open access repository and shared here. Please cite its original reference Running et al. (2004).

For reproducing results of Stocker et al. (2018b) regarding site-scale evaluations, run for example the scripts `plot_bias_all.R` and `plot_bias_problem.R`, available from (Github)[https://github.com/stineb/soilm_global] or (Zenodo)[https://zenodo.org/record/1286966#.W6TFipMzbUI], using CSV files provided here (see comments in scripts). For more insight, including analysis of global simulation outputs, see RMarkdown file `si_soilm_global.Rmd`. This renders the supplementary information PDF document provided along with Stocker et al. (2018b), which is available also on (RPubs)[http://rpubs.com/stineb/si_soilm_global2].

The datasets are prepared for open access distribution by script `prepare_data_openaccess.R ` on (Github)[https://github.com/stineb/soilm_global] or (Zenodo)[https://zenodo.org/record/1286966#.W6TFipMzbUI].

## Data description

### Site-level data

Data is provided as CSV files:

- `gpp_daily_fluxnet_stocker18natgeo.csv`: Daily data for full time series (not including MODIS GPP)
- `gpp_8daily_fluxnet_stocker18natgeo.csv`: Data aggregated to 8-day periods corresponding to MODIS dates (including MODIS GPP)
- `gpp_alg_daily_fluxnet_stocker18natgeo.csv`: Data filtered to periods with substantial soil moisture effects ("fLUE droughts" following Stocker et al. (2018a))
- `gpp_alg_8daily_fluxnet_stocker18natgeo.csv`: Data aggregated to 8-day periods and filtered to periods with substantial soil moisture effects.

Each column is a variable with the following name and units:

- `site_id`: FLUXNET site ID 
- `date`: Date of measurement, units: YYYY-MM-DD
- `gpp_obs`: Observed GPP from FLUXNET 2015 (see Stocker et al. (2018b), Methods, Observational Data), units: g C m-2 d-1
- `gpp_pmodel`: Simulated GPP from the P-model (see Stocker et al. (2018b), Methods, RS models), units: g C m-2 d-1
- `aet_splash`: Simulated actual evapotranspiration from the SPLASH model (Davis et al., 2017), units: mm d-1
- `pet_splash`: Simulated potential evapotranspiration from the SPLASH model (Davis et al., 2017), units: mm d-1
- `ppfd`: Photosynthetic photon flux density, based on shortwave incoming radiation data from FLUXNET 2015 (SW_IN_F), converted to PPFD as ppfd = SW_IN_F * kfFEC * 1.0e-6, with kfFEC = 2.04 micro-mol/J. Units: mol/m2/d
- `vpd`: Vapour pressure deficit, units: Pa
- `fpar_modis`: Fraction of absorbed photosynthetically active radiation by MODIS FPAR  (MCD15A3H Version 6, 500 m, 4-day) extracted using Google Earth Engine and the ‘gee_subset’ library (Hufkens, 2017), unitless.
- `soilm_mean`: Mean soil moisture used for analyses. This is based on observed soil moisture from FLUXNET 2015, mean across multiple depths if available, and normalised to between its minimum and maximum daily value for each site respectively, unitless.

### Global P-model simulation outputs

GPP and soil moisture output is provided as NetCDF files for simulations s0, s1a, s1b, and s1c (see Stocker et al. (2018b)). All meta information is provided therein. Files for simulation s1a are names as follows (for outputs from other simulations replace s1a with other simulation name).  

- `s1a_fapar3g_v2_global.d.gpp.nc`: Daily GPP from simulation s1a.
- `s1a_fapar3g_v2_global.d.wcont.nc`: Daily soil moisture from simulation s1a.


## References

Davis, T. W. et al. Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture. Geoscientific Model Development 10, 689–708 (2017).
Hufkens, K. khufkens/gee_subset: Google Earth Engine subset script & library. (2017). doi:10.5281/zenodo.833789Running, S. W. et al. A Continuous Satellite-Derived Measure of Global Terrestrial Primary Production. Bioscience 54, 547–560 (2004).
Stocker, B. et al., Quantifying soil moisture impacts on light use efficiency across biomes, New Phytologist, doi: 10.1111/nph.15123 (2018a).
Stocker, B. et al., Satellite monitoring underestimates the impact of drought on terrestrial primary productivity, Nature Geoscience (2018b). XXX STILL IN REVIEW XXX
Wang, H. et al. Towards a universal model for carbon dioxide uptake by plants. Nat Plants 3, 734–741 (2017).
