import  os
import  pylab            as pl
import  numpy            as np
import  redrock.io

from    redrock.external import desi
from    astropy.table    import Table, join


## export RR_TEMPLATE_DIR=$(pwd)/redrock-templates
os.environ["RR_TEMPLATE_DIR"] = "/Users/M.J.Wilson/work/cmbcross/mjw/desihub/redrock/py/redrock/templates"  

if __name__ == "__main__":
    '''
    This page documents how to run the `redrock.PlotSpec` viewer to get plots like:

    ![redrock.PlotSpec screenshot](rrplotspec.png)

    ## Installation
    Install the redrock + desispec code into an isolated conda environment:
    
    conda create -n rrdesi python=3 numpy scipy astropy numba ipython h5py matplotlib pyyaml
    
    pip install speclite

    source activate rrdesi

    git clone https://github.com/desihub/redrock-templates
    git clone https://github.com/desihub/redrock
    git clone https://github.com/desihub/desiutil
    git clone https://github.com/desihub/desispec

    for repo in desiutil desispec redrock; do
        cd $repo
        python setup.py install
        cd ..
        done

    ## Getting example data

    Available from NERSC at

    /project/projectdirs/desi/datachallenge/dc17a-twopct/dc17a-lite.tar.gz
    '''
    
    #- read redshift results and truth table
    specfile           = 'spectro/redux/dc17a2/spectra-64/172/17242/spectra-64-17242.fits'
    rr_h5file          = 'spectro/redux/dc17a2/spectra-64/172/17242/rr-64-17242.h5'
    zbestfile          = 'spectro/redux/dc17a2/spectra-64/172/17242/zbest-64-17242.fits'
    truthfile          = 'targets/truth-lite.fits'

    templates          = redrock.io.read_templates()

    targets, targetids = desi.read_spectra(specfile)

    zscan, zfit        = redrock.io.read_zscan(rr_h5file)
    
    zbest              = Table.read(zbestfile)  ## astropy tables
    truth              = Table.read(truthfile)

    ## p               = redrock.PlotSpec(targets, templates, zscan, zfit, truth=truth)

    #- Include truth info in plots; first fix some datamodel inconsistencies
    truth['targetid']  = truth['TARGETID']
    truth['ztrue']     = truth['TRUEZ']

    ## p               = redrock.PlotSpec(targets, templates, zscan, zfit, truth=truth)

    #- Filter out bad QSO targets
    ztruth                     = join(truth, zbest, keys='TARGETID')   ## join zbest and truth tables based on target id. 
    ztruth['TEMPLATETYPE'][:]  = np.char.strip(ztruth['TEMPLATETYPE'])

    dv                         = 3e5*(ztruth['Z'] - ztruth['TRUEZ']) / (1 + ztruth['TRUEZ'])

    isqso                      = ztruth['TEMPLATETYPE'] == 'QSO'

    bad                        = (np.abs(dv) > 1000) & (ztruth['ZWARN'] == 0) & (isqso)

    bad_targetids              = ztruth['TARGETID'][bad]
    bad_targets                = list()

    for t in targets:
        if t.id in bad_targetids:
            bad_targets.append(t)

    ## p                       = redrock.PlotSpec(bad_targets, templates, zscan, zfit, truth=truth)

    ## pl.savefig('rrplotspec2.png')
