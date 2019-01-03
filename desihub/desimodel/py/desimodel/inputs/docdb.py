'''
Utility functions for working with DocDB files
'''
import os

import numpy as np

from ..io import datadir

def _xls_col2int(col):
    '''
    Convert column string name to index, starting at 0
    
    e.g. A -> 0, B -> 1, ... Z -> 25, AA -> 26, AB -> 27
    '''
    abc = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    index = 0
    for i, x in enumerate(col.upper()[-1::-1]):
        index += (abc.index(x) + 1) * len(abc)**i

    return index-1

def xls_read_row(filename, sheetname, rownum, firstcol, lastcol, dtype=None):
    '''
    Read Excel file row from firstcol to lastcol
    
    Args:
        filename (str): Excel filename
        sheetname (str): sheet name within the filename
        rownum (int): 1-indexed row to read
        firstcol (str): Excel-style column name, e.g. 'A', 'B', or 'AC'
        lastcol (str): last column to include

    Options:
        dtype: convert output to this numpy dtype

    Returns numpy array of data
    
    Example:
        B5:D5 -> rownum=5, firstcol='B', lastcol='D' -> length 3 array
    '''
    import xlrd
    icol = _xls_col2int(firstcol)
    jcol = _xls_col2int(lastcol)
    with xlrd.open_workbook(filename) as wb:
        s = wb.sheet_by_name(sheetname)
        values = s.row_values(rownum-1, icol, jcol+1)
        return np.array(values, dtype=dtype)

def xls_read_col(filename, sheetname, column, firstrow, lastrow, dtype=None):
    '''
    Read Excel file column from firstrow to lastrow
    
    Args:
        filename (str): Excel filename
        sheetname (str): sheet name within the filename
        column (str): Excel-style column string, e.g. 'A', 'B', or 'AC'
        firstrow (int): 1-indexed first row to include
        lastrow (int): 1-indexed last row to include

    Options:
        dtype: convert output to this numpy dtype

    Returns numpy array of data
    
    Example:
        B5:B10 -> column='B', firstrow=5, lastrow=10 -> length 6 array
    '''
    import xlrd
    icol = _xls_col2int(column)
    with xlrd.open_workbook(filename) as wb:
        s = wb.sheet_by_name(sheetname)
        values = s.col_values(icol, firstrow-1, lastrow)
        return np.array(values, dtype=dtype)


def download(docnum, docver, filename, outdir=None, overwrite=False):
    '''
    Downloads and writes outdir/DESI-{docnum}v{docver}-{filename}

    Args:
        docnum: integer DocDB number
        docver: integer version number
        filename: string filename within that DocDB entry

    Options:
        outdir: output directory; default $DESIMODEL/data/inputs/docdb/
        overwrite: overwrite pre-existing file

    Returns:
        path to output file written

    Notes:

      * only supports python3
      * creates outdir if needed
      * prepends DESI-{docnum}v{docver} to {filename} even if filename
        already starts with that (in DocDB, some do and some don't...)      
    '''
    import urllib
    import requests
    from desiutil.log import get_logger
    log = get_logger()

    if outdir is None:
        outdir = os.path.join(datadir(), 'inputs', 'docdb')

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    outfile = 'DESI-{:04d}v{:d}-{}'.format(docnum, docver, filename)
    outfile = os.path.join(outdir, outfile)

    if os.path.exists(outfile):
        if overwrite:
            log.info('Redownloading and overwriting {}'.format(outfile))
        else:
            log.info('{} already exists; use overwrite=True to force redownload'.format(outfile))
            return outfile

    # e.g. https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=347;filename=DESI-347-v11%20Throughput%20Noise%20SNR%20Calcs.xlsx;version=11
    urlbase = 'https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile'
    url = '{}?docid={};version={};filename={}'.format(
        urlbase, docnum, docver, urllib.parse.quote(filename)
    )
    r = requests.get(url, auth=_auth())

    #- DocDB doesn't return the correct status codes for failure
    #- https://github.com/ericvaandering/DocDB/issues/11
    #- but just in case they do fix it in the future...
    if r.status_code != 200:
        raise IOError('Unable to download {}'.format(url))

    #- Work around what should have been a 404 Not Found
    failmsg = bytes('{} does not exist.'.format(filename), encoding='ascii')
    if b'There was a problem.' in r.content and failmsg in r.content:
        raise IOError('Unable to download {}'.format(url))

    with open(outfile, 'wb') as fx:
        fx.write(r.content)

    log.info('Wrote {}'.format(outfile))

    return outfile

#- lightly modified from desispec.download._auth;
#- consider refactoring into desiutil
def _auth(machine='desi.lbl.gov'):
    """Get authentication credentials.
    """
    from netrc import netrc
    from requests.auth import HTTPDigestAuth
    n = netrc()
    try:
        u,foo,p = n.authenticators(machine)
    except:
        raise ValueError('Unable to get user/pass from $HOME/.netrc for {}'.format(machine))

    return HTTPDigestAuth(u,p)    
