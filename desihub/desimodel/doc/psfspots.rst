=========
PSF Spots
=========

Overview
========

Spectrograph PSF files are in data/specpsf/desi-psf-\*.fits .  These
are in the `Specter "SpotGrid" format`_.

.. _`Specter "SpotGrid" format`: https://github.com/sbailey/specter/blob/dev/doc/datamodel/psf.md

Converting Zemax spots to Specter spots
=======================================

DocDB DESI-0334v1 contains spots from zemax + diffraction + CCD effects,
sampled on a grid of wavelength and slit position.  Here are the commands
to download and convert these to the files in etc/data/specpsf/.

YOU NORMALLY SHOULD NOT NEED TO RUN THESE.

They are for reference for when we need to generate new PSF files from a
new set of spots in DocDB::

    cd somewhere
    mkdir blue red nir

    DOCDB_USER=StephenBailey
    DOCDB_PASS=NotMyRealPassword

    cd blue
    SPOTFILE=DESI-0334-blue-images.zip
    wget --user $DOCDB_USER --password $DOCDB_PASS -O $SPOTFILE \
        "https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=334;filename=$SPOTFILE;version=1"
    unzip $SPOTFILE
    rm $SPOTFILE

    cd ../red
    SPOTFILE=DESI-0334-red-images.zip
    wget --user $DOCDB_USER --password $DOCDB_PASS -O $SPOTFILE \
        "https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=334;filename=$SPOTFILE;version=1"
    unzip $SPOTFILE
    rm $SPOTFILE

    cd ../nir
    SPOTFILE=DESI-0334-NIR-images-500.zip
    wget --user $DOCDB_USER --password $DOCDB_PASS -O $SPOTFILE \
        "https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=334;filename=$SPOTFILE;version=1"
    unzip $SPOTFILE
    rm $SPOTFILE

    cd ..
    python $DESIMODEL/bin/spots2psf.py blue/Blue*.fits --camera b \
        -o $DESIMODEL/data/specpsf/psf-b.fits
    python $DESIMODEL/bin/spots2psf.py red/Red*.fits --camera r \
        -o $DESIMODEL/data/specpsf/psf-r.fits
    python $DESIMODEL/bin/spots2psf.py nir/NIR*.fits --camera z \
        -o $DESIMODEL/data/specpsf/psf-z.fits

    rm -r blue nir red
