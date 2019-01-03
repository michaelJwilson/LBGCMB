export LEPHAREDIR="/Users/M.J.Wilson/work/cmbcross/mjw/dropouts"
export LEPHAREWORK="/Users/M.J.Wilson/work/cmbcross/mjw/dropouts/results/" 

## export PARAMS="/Users/M.J.Wilson/work/cmbcross/mjw/dropouts/para/original_zphot.para"
export PARAMS="/Users/M.J.Wilson/work/cmbcross/mjw/dropouts/para/zphot.para"
## export PARAMS="/Users/M.J.Wilson/work/cmbcross/mjw/dropouts/para/mjw_zphot.para"

## -- Configuration file
## emacs -nw /Users/M.J.Wilson/work/cmbcross/mjw/dropouts/zphot.para

### --- Spectra ---
## Binary Library for Stars
$LEPHAREDIR/source/sedtolib -c $PARAMS -t S

## Binary Library for QSO
$LEPHAREDIR/source/sedtolib -c $PARAMS -t Q

## Binary Library for Galaxies
## Check the GALAXY models used with: more $LEPHAREWORK/lib_bin/LIB_XXX.doc
$LEPHAREDIR/source/sedtolib -c $PARAMS -t G

### --- FILTERS ---
## Create a unique file with all transmission curves to be used
$LEPHAREDIR/source/filter -c $PARAMS

## Info about filters used
$LEPHAREDIR/source/filter_info -f HDF.filt -c $PARAMS

## Plot the transmission curves with a kshell including Supermongo macro
## $LEPHAREDIR/tools/filterplot $LEPHAREWORK/filt/HDF.filt

## Atmospheric and galactic extinction through the filters
## Cardelli law is hardcoded in the programs and is the default law for the galactic extinction; overwritten with -g; -e for atmospheric extinction.
$LEPHAREDIR/source/filter_extinc -f HDF.filt -o exti_test.out -g calzetti.dat -e extinc_eso.dat

### --- THEORETICAL MAGNITUDES AND K-CORRECTION COMPUTATION ---
## Mags for stars
$LEPHAREDIR/source/mag_star -c $PARAMS -MAGTYPE AB

## Mags for Gals                 
## Compute galaxy magnitudes and apply extinction for models ** between ** 4 to 8 with E(B-V)values read from the configuration file.                      
$LEPHAREDIR/source/mag_gal -c $PARAMS -t G

## Mags for QSOs; No extinction. 
## $LEPHAREDIR/source/mag_gal  -c $PARAMS -t Q -EB_V 0.0 -EM_LINES NO -EXTINC_LAW NONE

## Mags for Gals with fixed z-formation (if Z_FORM specified)
## $LEPHAREDIR/source/mag_zform -c $PARAMS

### --- PHOTO-Zs COMPUTATION ---
## Run the zphot code
## $LEPHAREDIR/source/zphota -c $PARAMS

## --- Supermongo plotting ---
## sm> macro read "{$LEPHAREDIR}/tools/spec.sm" zsp Name [1-131] [0-1]
## Name:  File name of the object to plot.
## [1-131] : SED of the best Gal/GAL-Far Infra Red/GAL-STOCHASTIC/QSO/Star to plot
## 1 : GAL
## 2 : GAL-FIR
## 4 : GAL-BCSTOCH
## 8 : QSO
## 16: STAR
## 32: GAL + GAL-FIR
## 64: GAL-BCSTOCH + GAL-FIR
## [1-131] Sum of any combinations of the above numbers
## [0-1]: 0 plot on screen, 1 plot in ps file
