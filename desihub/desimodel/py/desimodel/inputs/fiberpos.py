'''
Utilities for updating positioner to fiber mapping
'''
import os
import shutil

import numpy as np
from astropy.table import Table, vstack

from . import docdb
from ..io import datadir

def update(testdir=None, seed=2):
    '''
    Update positioner to fiber number mapping from DocDB

    Options:
        testdir: if not None, write files here instead of
            $DESIMODEL/data/footprint/fiberpos*
        seed: integer random number seed for randomization within a cartridge

    Writes testdir/fiberpos* or $DESIMODEL/data/focalplane/fiberpos*
    '''
    from desiutil.log import get_logger
    log = get_logger()

    #- Download input files from DocDB
    cassette_file = docdb.download(2721, 2, 'cassette_order.txt')
    xls_fp_layout = docdb.download(530, 11, 'DESI-0530-v11 (Focal Plane Layout).xlsx')
    platescale_file = docdb.download(329, 15, 'Echo22Platescale.txt')

    #- Pick filenames in output directory
    if testdir is None:
        outdir = os.path.join(datadir(), 'focalplane')
    else:
        outdir = testdir

    if not os.path.isdir(outdir):
        raise ValueError("Missing directory {}".format(testdir))

    #- copy platescale file
    outpsfile = os.path.join(outdir, 'platescale.txt')
    shutil.copy(platescale_file, outpsfile)
    log.info('Wrote {}'.format(outpsfile))

    #- Random but reproducible
    np.random.seed(seed)

    #- DESI-0530 file name (fn) and sheet name (sn) shortcuts
    fn = xls_fp_layout
    sn = 'PositionerAndFiducialLocations'

    #- Sanity check that columns are still in the same place
    rowmin, rowmax = 48, 590
    headers = docdb.xls_read_row(fn, sn, rowmin-1, 'B', 'S')

    assert headers[0] == 'device_location_id'
    assert headers[1] == 'device_type'
    assert headers[2] == 'X'
    assert headers[3] == 'Y'
    assert headers[4] == 'Z'
    assert headers[8] == 'cassetteID'
    assert headers[15] == 'Q'
    assert headers[17] == 'S'

    #- Read Excel table with device locations
    posloc = Table()
    posloc['DEVICE'] = docdb.xls_read_col(fn, sn, 'B', rowmin, rowmax, dtype=int)
    posloc['DEVICE_TYPE'] = docdb.xls_read_col(fn, sn, 'C', rowmin, rowmax, dtype=str)
    posloc['X'] = docdb.xls_read_col(fn, sn, 'D', rowmin, rowmax, dtype=float)
    posloc['Y'] = docdb.xls_read_col(fn, sn, 'E', rowmin, rowmax, dtype=float)
    posloc['Z'] = docdb.xls_read_col(fn, sn, 'F', rowmin, rowmax, dtype=float)
    posloc['Q'] = docdb.xls_read_col(fn, sn, 'Q', rowmin, rowmax, dtype=float)
    posloc['S'] = docdb.xls_read_col(fn, sn, 'S', rowmin, rowmax, dtype=float)

    #- Cassette N/A -> -1, and parse string -> float -> int
    c = docdb.xls_read_col(fn, sn, 'J', rowmin, rowmax)
    not_spectro_fiber = (c == 'N/A')
    c[not_spectro_fiber] = '-1'
    posloc['CASSETTE'] = np.array(c, dtype=float).astype(int)

    #- Sanity check on values
    ndevice = len(posloc)
    assert ndevice == 543  #- 543 holes have been drilled
    assert len(np.unique(posloc['DEVICE'])) == len(posloc['DEVICE'])
    assert set(posloc['DEVICE_TYPE']) == set(['POS', 'FIF', 'GIF', 'NON', 'OPT', 'ETC'])
    assert 0 < np.min(posloc['X']) and np.max(posloc['X']) < 410
    assert 0 <= np.min(posloc['Q']) and np.max(posloc['Q']) < 36.0
    assert 0 <= np.min(posloc['S']) and np.max(posloc['S']) < 412.3
    assert np.all(posloc['S']**2 > posloc['X']**2 + posloc['Y']**2 + posloc['Z']**2)
    assert np.min(posloc['CASSETTE']) == -1
    assert np.max(posloc['CASSETTE']) == 11
    assert set(posloc['DEVICE_TYPE'][posloc['CASSETTE']==11]) == set(['ETC', 'OPT'])
    assert set(posloc['DEVICE_TYPE'][posloc['CASSETTE']==-1]) == set(['FIF', 'GIF', 'NON'])
    assert 0 not in posloc['CASSETTE']

    #- Read mapping of cassettes on focal plane to fibers on slithead
    colnames = ['fibermin', 'fibermax', 'sp0', 'sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8', 'sp9']
    cassettes = Table.read(cassette_file, format='ascii', names=colnames)

    #- Randomize fibers within a cassette
    petals = list()
    for p in range(10):
        fiberpos = posloc.copy(copy_data=True)
        fiberpos['FIBER'] = -1
        fiberpos['PETAL'] = p
        fiberpos['SLIT'] = p
        fiberpos['SPECTRO'] = p
        iipos = (fiberpos['DEVICE_TYPE'] == 'POS')
        ### fiberpos['device'] += p*len(fiberpos)
        for c in range(1,11):
            ii = (cassettes['sp'+str(p)] == c)
            assert np.count_nonzero(ii) == 1
            fibermin = p*500 + cassettes['fibermin'][ii][0]
            fibermax = p*500 + cassettes['fibermax'][ii][0]

            jj = iipos & (fiberpos['CASSETTE'] == c)
            assert np.count_nonzero(jj) == 50
            fiber = list(range(fibermin, fibermax+1))
            np.random.shuffle(fiber)
            fiberpos['FIBER'][jj] = fiber

        #- Additional columns
        fiberpos['SLITBLOCK'] = (fiberpos['FIBER'] % 500) // 25
        fiberpos['BLOCKFIBER'] = (fiberpos['FIBER'] % 500) % 25
        fiberpos['LOCATION'] = p*1000 + fiberpos['DEVICE']

        #- Petal 0 is at the "bottom"; See DESI-0530
        phi = np.radians((7*36 + 36*p)%360)
        x = np.cos(phi)*fiberpos['X'] - np.sin(phi)*fiberpos['Y']
        y = np.sin(phi)*fiberpos['X'] + np.cos(phi)*fiberpos['Y']
        fiberpos['X'] = x
        fiberpos['Y'] = y

        petals.append(fiberpos)

    fiberpos = vstack(petals)
    fiberpos.sort('FIBER')
    POS = (fiberpos['DEVICE_TYPE'] == 'POS')

    #- devices that don't go to spectrographs don't have slitblock, blockfiber
    fiberpos['SLITBLOCK'][~POS] = -1
    fiberpos['BLOCKFIBER'][~POS] = -1

    #- More sanity checks before writing output
    fp = fiberpos[POS]
    assert len(fp) == 5000
    assert len(np.unique(fp['FIBER'])) == 5000
    assert min(fp['FIBER']) == 0
    assert max(fp['FIBER']) == 4999
    assert len(set(fp['SPECTRO'])) == 10
    assert min(fp['SPECTRO']) == 0
    assert max(fp['SPECTRO']) == 9
    assert len(np.unique(fiberpos['DEVICE'])) == ndevice
    assert len(np.unique(fiberpos['LOCATION'])) == len(fiberpos)

    #- Drop some columns we don't need
    fiberpos.remove_column('CASSETTE')

    #- Update i8 -> i4 for integer columns
    for colname in ['FIBER', 'DEVICE', 'SPECTRO', 'PETAL', 'SLIT']:
        fiberpos.replace_column(colname, fiberpos[colname].astype('i4'))

    #- Reorder columns
    assert set(fiberpos.colnames) == set('DEVICE DEVICE_TYPE X Y Z Q S FIBER PETAL SLIT SPECTRO SLITBLOCK BLOCKFIBER LOCATION'.split())
    colnames = 'PETAL DEVICE DEVICE_TYPE LOCATION FIBER X Y Z Q S  SPECTRO SLIT SLITBLOCK BLOCKFIBER'.split()
    fiberpos = fiberpos[colnames]
    assert fiberpos.colnames == colnames

    #- Set units and descriptions; see DESI-2724
    fiberpos['X'].unit = 'mm'
    fiberpos['Y'].unit = 'mm'
    fiberpos['Z'].unit = 'mm'
    fiberpos['Q'].unit = 'deg'
    fiberpos['S'].unit = 'mm'
    fiberpos['X'].description = 'focal surface location [mm]'
    fiberpos['Y'].description = 'focal surface location [mm]'
    fiberpos['Z'].description = 'focal surface location [mm]'
    fiberpos['Q'].description = 'azimuthal angle on focal surface [deg]'
    fiberpos['S'].description = 'radial distance along focal surface [mm]'
    fiberpos['FIBER'].description = 'fiber number [0-4999]'
    fiberpos['DEVICE'].description = 'focal plane device_loc number [0-542]'
    fiberpos['SPECTRO'].description = 'spectrograph number [0-9]'
    fiberpos['PETAL'].description = 'focal plane petal_loc number [0-9]'
    fiberpos['SLIT'].description = 'spectrograph slit number [0-9]'
    fiberpos['SLITBLOCK'].description = 'id of the slitblock on the slit [0-19]'
    fiberpos['BLOCKFIBER'].description = 'id of the fiber on the slitblock [0-24]'
    fiberpos['LOCATION'].description = 'global location id across entire focal plane [0-9543]; has gaps in sequence'
    fiberpos.meta['comments'] = [
        "Coordinates at zenith: +x = East = +RA; +y = South = -dec",
        "PETAL and DEVICE refer to locations, not hardware serial numbers",
        "Differences from DESI-2724 naming:",
        '  - Drops "_ID" from column names',
        '  - Drops "_LOC" from "DEVICE_LOC" and "PETAL_LOC"',
        "  - SLITBLOCK as int [0-19] instead of string [B0-B19]",
        "  - BLOCKFIBER as int [0-24] instead of string [F0-F24]",
        "Convenience columns:",
        "  - FIBER = PETAL*500 + SLITBLOCK*25 + BLOCKFIBER",
        "  - LOCATION = PETAL*1000 + DEVICE",
    ]

    ecsvout = os.path.join(outdir, 'fiberpos.ecsv')
    textout = os.path.join(outdir, 'fiberpos.txt')
    fitsout = os.path.join(outdir, 'fiberpos.fits')
    pngout  = os.path.join(outdir, 'fiberpos.png')

    #- Write old text format with just fiber, device, spectro, x, y, z
    write_text_fiberpos(textout, fiberpos[POS])
    log.info('Wrote {}'.format(textout))

    #- Write all columns but only for positioners with fibers
    fiberpos[POS].write(ecsvout, format='ascii.ecsv')
    log.info('Wrote {}'.format(ecsvout))

    fiberpos[POS].write(fitsout, format='fits', overwrite=True)
    log.info('Wrote {}'.format(fitsout))

    #- Write all columns and all rows, including
    #- fiducials (device_type='FIF') and sky monitor (device_type='ETC')
    fiberpos.sort('LOCATION')
    fitsallout = fitsout.replace('.fits', '-all.fits')
    ecsvallout = textout.replace('.txt', '-all.ecsv')
    fiberpos.write(fitsallout, format='fits', overwrite=True)
    fiberpos.write(ecsvallout, format='ascii.ecsv')
    log.info('Wrote {}'.format(fitsallout))
    log.info('Wrote {}'.format(ecsvallout))

    #- Visualize mapping
    POS = (fiberpos['DEVICE_TYPE'] == 'POS')
    FIF = (fiberpos['DEVICE_TYPE'] == 'FIF')
    ETC = (fiberpos['DEVICE_TYPE'] == 'ETC')
    import pylab as P
    P.jet()     #- With apologies to viridis
    P.figure(figsize=(7,7))
    P.scatter(fiberpos['X'][POS], fiberpos['Y'][POS], c=fiberpos['FIBER'][POS]%500, edgecolor='none', s=20)
    # P.scatter(fiberpos['x'][FIF], fiberpos['y'][FIF], s=5, color='k')
    # P.plot(fiberpos['x'][ETC], fiberpos['y'][ETC], 'kx', ms=3)
    P.grid(alpha=0.2, color='k')
    P.xlim(-420,420)
    P.ylim(-420,420)
    P.xlabel('x [mm]')
    P.ylabel('y [mm]')
    P.title('Focal plane color coded by fiber location on slithead')
    P.savefig(pngout, dpi=80)
    log.info('Wrote {}'.format(pngout))

def write_text_fiberpos(filename, fiberpos):
    '''
    Writes a fiberpos table to filename, maintaining backwards compatibility
    with the original fiberpos.txt format

    Args:
        filename: output file name string
        fiberpos: astropy Table of fiber positions
    '''

    #- Write the old text file format for backwards compatibility
    fxlines = [
        "#- THIS FILE IS PROVIDED FOR BACKWARDS COMPATIBILITY",
        "#- Please use fiberpos-all.[ecsv,fits] for additional columns",
        '#- and non-spectrofiber device hole locations.',
        '#-'
        "#- Fiber to focal plane device hole mapping; x,y,z in mm on focal plane",
        "#- See doc/fiberpos.rst and DESI-0530 for more details.",
        "#- Coordinates at zenith: +x = East = +RA; +y = South = -dec",
        "",
        "#- fiber=at spectrograph; fpdevice=numbering on focal plane",
        "",
        '#- fiber location spectro  x  y  z']
    for row in fiberpos:
        fxlines.append("{:4d}  {:4d}  {:2d}  {:12.6f}  {:12.6f}  {:12.6f}".format(
            row['FIBER'], row['LOCATION'], row['SPECTRO'],
            row['X'], row['Y'], row['Z'],
        ))

    with open(filename, 'w') as fx:
        fx.write('\n'.join(fxlines)+'\n')
