# Data from article Stocker et al. (2018b) *Nature Geosci.*

The datasets provided here include:

- Site-level GPP model results from the P-model (Wang et al., 2017)
- Model outputs from global simulations with the P-model (Wang et al., 2017) as implemented for the study by Stocker et al. (2018b)

This data may be used to partly reproduce results presented in Stocker et al. (2018b) *Nature Geosci.*. "Partly" because we used data for our analysis that was not open access but was confidentially shared with us. This includes remote sensing-based GPP estimates from the BESS and VPM models. Other open access data that was used for the analysis may not be distributed under this DOI. This includes FLUXNET 2015 data and MODIS data.

For reproducing results of Stocker et al. (2018b) regarding site-scale evaluations, run for example the scripts `plot_bias_all.R` and `plot_bias_problem.R`, available from (Github)[https://github.com/stineb/soilm_global] or (Zenodo)[https://zenodo.org/record/1286966#.W6TFipMzbUI], using CSV files provided here (see comments in scripts). For more insight, including analysis of global simulation outputs, see RMarkdown file `si_soilm_global.Rmd`. This renders the supplementary information PDF document provided along with Stocker et al. (2018b), which is available also on (RPubs)[http://rpubs.com/stineb/si_soilm_global2].

The present datasets are prepared by script `prepare_data_openaccess.R ` on (Github)[https://github.com/stineb/soilm_global] or (Zenodo)[https://zenodo.org/record/1286966#.W6TFipMzbUI].

## Data description

### Site-level data

Data is provided as CSV files:

- `gpp_daily_fluxnet_stocker18natgeo.csv`: Daily data for full time series (not including MODIS GPP)
- `gpp_8daily_fluxnet_stocker18natgeo.csv`: Data aggregated to 8-day periods corresponding to MODIS dates (including MODIS GPP)
- `gpp_alg_daily_fluxnet_stocker18natgeo.csv`: Data filtered to periods with substantial soil moisture effects ("fLUE droughts" following Stocker et al. (2018a))
- `gpp_alg_8daily_fluxnet_stocker18natgeo.csv`: Data aggregated to 8-day periods and filtered to periods with substantial soil moisture effects.

Each column is a variable with the following name and units (not all variables are available in all files):

- `site_id`: FLUXNET site ID 
- `date`: Date of measurement, units: YYYY-MM-DD
- `gpp_pmodel` and `gpp_modis`: Simulated GPP from the P-model and MODIS (see Stocker et al. (2018b), Methods, RS models), units: g C m-2 d-1 (mean across 8 day periods in respective files)
- `aet_splash`: Simulated actual evapotranspiration from the SPLASH model (Davis et al., 2017), units: mm d-1
- `pet_splash`: Simulated potential evapotranspiration from the SPLASH model (Davis et al., 2017), units: mm d-1
- `soilm_splash`: Soil moisture simulated by the SPLASH model (Davis et al., 2017), normalised to vary between zero and one at the maximum water holding capacity, unitless.
- `flue`: fLUE estimate from Stocker et al. (2018). Estimates soil moisture stress on light use efficiency from flux data, unitless.
- `beta_a`, `beta_b`, and `beta_c`: Empirical soil moisture stress, used as multiplier to simulated GPP as described in Stocker et al. (2018b), unitless.

### Global P-model simulation outputs

GPP and soil moisture output is provided as NetCDF files for simulations s0, and s1b (see Stocker et al. (2018b)). All meta information is provided therein. Files for simulation s1b are names as follows (for outputs from other simulations replace s1b with other simulation name). The fraction of each gridcell covered by land (not open water or ice) is given by separate file `s1b_fapar3g_v2_global.fland.nc`.

- `s1b_fapar3g_v2_global.d.gpp.nc`: Daily GPP from simulation s1b.
- `s1b_fapar3g_v2_global.d.wcont.nc`: Daily soil moisture from simulation s1b.

Due to limited total file size allowed for uploads to Zenodo, only outputs from s1b are provided here. Other outputs may be obtained upon request addressed to benjamin.stocker@gmail.com. 

## References

Davis, T. W. et al. Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture. Geoscientific Model Development 10, 689–708 (2017).
Hufkens, K. khufkens/gee_subset: Google Earth Engine subset script & library. (2017). doi:10.5281/zenodo.833789Running, S. W. et al. A Continuous Satellite-Derived Measure of Global Terrestrial Primary Production. Bioscience 54, 547–560 (2004).
Stocker, B. et al., Quantifying soil moisture impacts on light use efficiency across biomes, New Phytologist, doi: 10.1111/nph.15123 (2018a).
Stocker, B. et al., Satellite monitoring underestimates the impact of drought on terrestrial primary productivity, Nature Geoscience (2018b). XXX STILL IN REVIEW XXX
Wang, H. et al. Towards a universal model for carbon dioxide uptake by plants. Nat Plants 3, 734–741 (2017).
