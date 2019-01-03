=========
desimodel
=========

Introduction
------------

This product contains information about DESI hardware designs in machine
readable formats to be used by simulations.  It is intended to be correct
but non-authoritative.  If in question, the DocDB version of a design
parameter is correct.  This is intended to reflect the same information
while being more conveniently organized and formatted for simulations.

**PLEASE KEEP THIS FILE IN SYNC WITH THE EQUIVALENT FILE IN SVN.**

Desimodel Code
--------------

The svn product described below contains only the *data* associated with
desimodel. The code is in Github: https://github.com/desihub/desimodel.

Desimodel Data
--------------

Getting the Data
~~~~~~~~~~~~~~~~

The data that accompanies the desimodel code is not stored with the code.
Due to its size, it is kept in the DESI svn repository.  The public, read-only
URL for svn access is https://desi.lbl.gov/svn/code/desimodel, with the usual
trunk/, tags/ and branches/ directories.

Once you have installed this package, using either pip or desiInstall, there
are two ways to install the accompanying data.  **For most every case, you
should install the tag in svn that corresponds to the same tag in git.**

There are two methods to install the data, "by hand" and "scripted."

For "by hand" installs:

1. Define the environment variable ``DESIMODEL``::

    export DESIMODEL=/Users/desicollaborator/Data/desimodel/0.4.2

   Note how the tag name is included.

2. Create the directory and switch to it::

    mkdir -p $DESIMODEL
    cd $DESIMODEL

3. Export::

    svn export https://desi.lbl.gov/svn/code/desimodel/tags/0.4.2/data

   Note how the tag name is the *same* as in the ``DESIMODEL`` variable.

4. You may now want to add ``DESIMODEL`` to your shell startup scripts.

For "scripted" installs:

* Installing this package will create the command-line script
  ``install_desimodel_data``.  It should appear in your ``PATH`` assuming
  a successful install.  ``install_desimodel_data --help`` will show you
  how to use this script.  Basically it is just a wrapper on the "by hand"
  method described above.
* You can also call the function ``desimodel.install.install()`` from
  inside other Python code.

Regardless of which method you choose, you should set the ``DESIMODEL``
environment variable to point to the directory containing the data/
directory.  The only real difference among all these methods is exactly
*when* you define the ``DESIMODEL`` variable.

Data Files
~~~~~~~~~~

Files in data/inputs/ are copied from DocDB and are not intended to be used
directly.  These are parsed and reformatted to produce files in other data/
directories which are for use.

data/desi.yaml
    Basic scalar parameters, organized in a nested tree.

data/focalplane/
    Informaton about positioner locations and platescale.

data/specpsf/
    Spectrograph point-spread-function (PSF) for specter_
    CCD pixel-level simulations.

data/throughput/
    Throughput *versus* wavelength (also contained in specpsf).

data/footprint/
    DESI footprint with RA, Dec of tiles.

data/spectra/
    Example benchmark spectra.

.. _specter: https://github.com/desihub/specter

Branches
~~~~~~~~

There are a couple of permanent branches that were used for testing
alternative designs.  These will never be merged into trunk but we
will keep them around for the record:

altccd
    500 micron *versus* 250 micron thick CCDs.

bb
    Recreating assumptions used during early BigBOSS projections.

In addition to these historical branches, there is a permanent 'testing' branch
that contains smaller versions of the desimodel files.  This branch is
intended for use in desimodel unit tests.

Tagging
-------

If *either* the data *or* the code changes, a new tag should be created in
both git and svn.

Full Documentation
------------------

Please visit `desimodel on Read the Docs`_

.. image:: http://readthedocs.org/projects/desimodel/badge/?version=latest
    :target: http://desimodel.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. _`desimodel on Read the Docs`: http://desimodel.readthedocs.io/en/latest/

Travis Build Status
-------------------

.. image:: https://img.shields.io/travis/desihub/desimodel.svg
    :target: https://travis-ci.org/desihub/desimodel
    :alt: Travis Build Status


Test Coverage Status
--------------------

.. image:: https://coveralls.io/repos/desihub/desimodel/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/desihub/desimodel?branch=master
    :alt: Test Coverage Status

License
-------

desimodel is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
