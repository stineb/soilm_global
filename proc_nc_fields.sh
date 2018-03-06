#!/bin/bash

get_anom () {

	##-----------------------------
	## argument 1: base file name
	##-----------------------------
	echo "----------------------------"
	echo "processing $1 ..."
	echo "----------------------------"

	## detrend and remove mean	
	echo "detrending $1_MON.nc..."
	cdo detrend -selyear,1982/2016 -selname,gpp $1_MON.nc $1_DETR_MON_0mean.nc

	## get daily anomalies, warning: this doesn't remove the linear trend
	echo "get anomalies..."
	cdo -ymonsub $1_DETR_MON_0mean.nc -ymonmean $1_DETR_MON_0mean.nc $1_MON_ANOM.nc

	rm $1_DETR_MON_0mean.nc

	return 0

}


proc_30y () {
	##-----------------------------
	## argument 1: base file name
	##-----------------------------
	echo "----------------------------"
	echo "processing $1 ..."
	echo "----------------------------"

	## detrend at each gridcell
	cdo detrend -selyear,1982/2016 -selname,gpp $1_ANN.nc $1_DETR.nc
	cdo detrend -selyear,2001/2016 -selname,gpp $1_ANN.nc $1_DETR20XX.nc

	## get variance of annual GPP at each pixel
	cdo timvar $1_DETR.nc $1_VAR.nc
	cdo timvar $1_DETR20XX.nc $1_VAR20XX.nc

	## get mean field
	cdo timmean -selyear,1982/2016 -selname,gpp $1_ANN.nc $1_MEAN.nc
	cdo timmean -selyear,2001/2016 -selname,gpp $1_ANN.nc $1_MEAN20XX.nc

	## get relative variance field
	cdo div $1_VAR.nc $1_MEAN.nc $1_RELVAR.nc
	cdo div $1_VAR20XX.nc $1_MEAN20XX.nc $1_RELVAR20XX.nc

	## get global totals
	## GPP
	cdo gridarea $1_ANN.nc gridarea.nc
	cdo mulc,1 -seltimestep,1 $1_ANN.nc tmp.nc
	cdo div tmp.nc tmp.nc ones.nc
	cdo selname,gpp ones.nc mask.nc
	cdo mul mask.nc gridarea.nc gridarea_masked.nc
	cdo mul gridarea_masked.nc $1_ANN.nc tmp2.nc
	cdo fldsum tmp2.nc tmp3.nc
	cdo mulc,1e-15 tmp3.nc $1_GLOB.nc

	## additional file
	cdo selyear,1982/2016 $1_GLOB.nc $1_GLOB30yr.nc
	cdo selyear,2001/2016 $1_GLOB.nc $1_GLOB20XX.nc

	## detrend
	cdo detrend -selyear,1982/2016 -selname,gpp $1_GLOB.nc $1_DETR_GLOB.nc
	cdo detrend -selyear,2001/2016 -selname,gpp $1_GLOB.nc $1_DETR_GLOB20XX.nc

	## variance of global time series
	cdo timvar $1_DETR_GLOB.nc $1_VAR_GLOB.nc
	cdo timvar $1_DETR_GLOB20XX.nc $1_VAR_GLOB20XX.nc

	## mean of global time series
	cdo timmean $1_GLOB30yr.nc $1_MEAN_GLOB.nc
	cdo timmean $1_GLOB20XX.nc $1_MEAN_GLOB20XX.nc

	## relative variance of global time series
	cdo div $1_VAR_GLOB.nc $1_MEAN_GLOB.nc $1_RELVAR_GLOB.nc
	cdo div $1_VAR_GLOB20XX.nc $1_MEAN_GLOB20XX.nc $1_RELVAR_GLOB20XX.nc

	## remove temporary files
	rm tmp.nc tmp2.nc tmp3.nc gridarea.nc gridarea_masked.nc mask.nc ones.nc

	return 0

}

proc_10y () {
	##-----------------------------
	## argument 1: base file name
	##-----------------------------
	echo "----------------------------"
	echo "processing $1 ..."
	echo "----------------------------"

	## detrend at each gridcell
	cdo detrend -selyear,2001/2011 -selname,gpp $1_ANN.nc $1_DETR20XX.nc

	## get variance of annual GPP at each pixel
	cdo timvar $1_DETR20XX.nc $1_VAR20XX.nc

	## get mean field
	cdo timmean -selyear,2001/2011 -selname,gpp $1_ANN.nc $1_MEAN20XX.nc

	## get relative variance field
	cdo div $1_VAR20XX.nc $1_MEAN20XX.nc $1_RELVAR20XX.nc

	## get global totals
	## GPP
	cdo gridarea $1_ANN.nc gridarea.nc
	cdo mulc,1 -seltimestep,1 $1_ANN.nc tmp.nc
	cdo div tmp.nc tmp.nc ones.nc
	cdo selname,gpp ones.nc mask.nc
	cdo mul mask.nc gridarea.nc gridarea_masked.nc
	cdo mul gridarea_masked.nc $1_ANN.nc tmp2.nc
	cdo fldsum tmp2.nc tmp3.nc
	cdo mulc,1e-15 tmp3.nc $1_GLOB.nc

	## additional file
	cdo selyear,2001/2011 $1_GLOB.nc $1_GLOB20XX.nc

	## detrend
	cdo detrend -selyear,2001/2011 -selname,gpp $1_GLOB.nc $1_DETR_GLOB20XX.nc

	## variance
	cdo timvar $1_DETR_GLOB20XX.nc $1_VAR_GLOB20XX.nc

	## mean of global time series
	cdo timmean $1_GLOB20XX.nc $1_MEAN_GLOB20XX.nc

	## relative variance of global time series
	cdo div $1_VAR_GLOB20XX.nc $1_MEAN_GLOB20XX.nc $1_RELVAR_GLOB20XX.nc
	
	## remove temporary files
	rm tmp.nc tmp2.nc tmp3.nc gridarea.nc gridarea_masked.nc mask.nc ones.nc

	return 0

}

here=`pwd`
myhome=~


##----------------------------------------------------
## P-model S0
##----------------------------------------------------
cd $myhome/data/pmodel_fortran_output/

## reduce to 30 years and get annual total
echo "get annual totals..."
cdo monsum  -selyear,1982/2016 s0_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s0_MON.nc
cdo yearsum -selyear,1982/2016 s0_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s0_ANN.nc

## process annual files
echo "process annual files..."
proc_30y gpp_pmodel_s0

## process daily files to get daily anomalies
ln -svf s0_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s0_DAY.nc
get_anom gpp_pmodel_s0

cd $here


##----------------------------------------------------
## P-model S1a, S1b, S1c
##----------------------------------------------------
cd $myhome/data/pmodel_fortran_output/

## reduce to 30 years and get annual total
echo "get annual totals..."
cdo yearsum -selyear,1982/2016 s1a_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1a_ANN.nc
cdo yearsum -selyear,1982/2016 s1b_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1b_ANN.nc
cdo yearsum -selyear,1982/2016 s1c_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1c_ANN.nc

cdo monsum  -selyear,1982/2016 s1a_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1a_MON.nc
cdo monsum  -selyear,1982/2016 s1b_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1b_MON.nc
cdo monsum  -selyear,1982/2016 s1c_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1c_MON.nc

## process annual files
echo "process annual files..."
proc_30y gpp_pmodel_s1a
proc_30y gpp_pmodel_s1b
proc_30y gpp_pmodel_s1c

echo "get ensemble means..."

## get mean annual GPP across S1a, S1b, S1c
cdo ensmean gpp_pmodel_s1a_MEAN.nc gpp_pmodel_s1b_MEAN.nc gpp_pmodel_s1c_MEAN.nc gpp_pmodel_s1_MEAN.nc 

## get mean variance S1a, S1b, S1c
cdo ensmean gpp_pmodel_s1a_VAR.nc gpp_pmodel_s1b_VAR.nc gpp_pmodel_s1c_VAR.nc gpp_pmodel_s1_VAR.nc

## get mean relative variance S1a, S1b, S1c
cdo ensmean gpp_pmodel_s1a_RELVAR.nc gpp_pmodel_s1b_RELVAR.nc gpp_pmodel_s1c_RELVAR.nc gpp_pmodel_s1_RELVAR.nc

## get daily soil moisture limitation factor from comparing GPP to s0
echo "get daily soil moisture limitation factor..."
cdo div s1a_fapar3g_v2_global.d.gpp.nc s0_fapar3g_v2_global.d.gpp.nc s1a_fapar3g_v2_global.d.soilmstress.nc
cdo div s1b_fapar3g_v2_global.d.gpp.nc s0_fapar3g_v2_global.d.gpp.nc s1b_fapar3g_v2_global.d.soilmstress.nc
cdo div s1c_fapar3g_v2_global.d.gpp.nc s0_fapar3g_v2_global.d.gpp.nc s1c_fapar3g_v2_global.d.soilmstress.nc

## process daily files to get daily anomalies
ln -svf s1a_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1a_DAY.nc
ln -svf s1b_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1b_DAY.nc
ln -svf s1c_fapar3g_v2_global.d.gpp.nc gpp_pmodel_s1c_DAY.nc

get_anom gpp_pmodel_s1a
get_anom gpp_pmodel_s1b
get_anom gpp_pmodel_s1c

cd $here


# ##----------------------------------------------------
# ## MTE
# ##----------------------------------------------------
# cd $myhome/data/gpp_mte/

# ## get a new time axis (file starts in december 1981 when it should be jan 1982)
# Rscript preproc_mte.R

# ## select years
# cdo selyear,1982/2011 gpp_mte_NICE.nc gpp_mte_SUB.nc

# ## multiply with days per month
# cdo muldpm gpp_mte_SUB.nc gpp_mte_DPM.nc

# ## multiply with seconds per day and convert from kg C to g C
# cdo mulc,86400000 gpp_mte_DPM.nc gpp_mte_SPM.nc

# ## get annual sums
# cdo yearsum gpp_mte_SPM.nc gpp_mte_ANN.nc

# ## process fully
# proc_30y gpp_mte

# cd $here


# ##----------------------------------------------------
# ## MTE-FLUXCOM
# ##----------------------------------------------------
# cd $myhome/data/gpp_mte/

# ## concatenae annual files (cdo mergetime doesn't work)
# Rscript preprocess_gpp_mte_fluxcom.R

# proc_30y gpp_mte_fluxcom

# cd $here

# ##----------------------------------------------------
# ## VPM
# ##----------------------------------------------------
# cd $myhome/data/gpp_vpm/

# # ## concatenae annual files (cdo mergetime doesn't work)
# # Rscript preprocess_vpm.R

# ## process fully
# proc_10y gpp_vpm

# cd $here


# ##----------------------------------------------------
# ## BESS
# ##----------------------------------------------------
# cd $myhome/data/gpp_bess/

# ## multiply with days per month
# cdo muldpm gpp_bess.nc gpp_bess_DPM.nc

# ## get annual sums
# cdo yearsum gpp_bess_DPM.nc gpp_bess_ANN.nc

# ## process fully
# proc_10y gpp_bess

# cd $here


# ##----------------------------------------------------
# ## MODIS
# ##----------------------------------------------------
# cd $myhome/data/gpp_modis/

# ## process fully
# proc_10y gpp_modis

# cd $here

















