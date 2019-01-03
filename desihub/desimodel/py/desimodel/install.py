# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.install
=================

Install data files not handled by pip install.
"""


def default_install_dir():
    """Return the default install directory.  Assumes this file lives in
    a 'site-packages' directory.

    Returns
    -------
    :class:`str`
        The path to the install directory.
    """
    from os.path import dirname
    return dirname(dirname(dirname(dirname(dirname(__file__)))))


def assert_svn_exists():
    """Assert svn command exists and raise an informative error if not"""

    from subprocess import check_output, CalledProcessError
    try:
        r = check_output(['svn', '--version'])
    except OSError as e:
        raise AssertionError("svn command is not executable. Install svn to use the install script. Original Error is %s", e.strerror)
    except CalledProcessError as e:
        raise AssertionError("The svn command (%s) on this system does not work. Output is %s" % (e.cmd, e.output))

def svn_export(desimodel_version=None):
    """Create a :command:`svn export` command suitable for downloading a
    particular desimodel version.

    Parameters
    ----------
    desimodel_version : :class:`str`, optional
        The version X.Y.Z to download, trunk, or something of the
        form branches/... Defaults to trunk.

    Returns
    -------
    :class:`list`
        A :command:`svn` command in list form, suitable for passing to
        :class:`subprocess.Popen`.
    """
    from . import __version__ as this_version
    if desimodel_version is None:
        export_version = 'trunk'
    elif desimodel_version is 'trunk' or 'branches/' in desimodel_version:
        export_version = desimodel_version
    else:
        export_version = 'tags/' + desimodel_version
    return ["svn", "export",
            ("https://desi.lbl.gov/svn/code/desimodel/" +
             "{0}/data").format(export_version)]


def install(desimodel=None, version=None):
    """Primary workhorse function.

    Parameters
    ----------
    desimodel : :class:`str`, optional
        Allows the install directory to be explicitly set.
    version : :class:`str`, optional
        Allows the desimodel version to be explicitly set.

    Raises
    ------
    :class:`RuntimeError`
        Standard error output from svn export command when status is non-zero.
    """
    from os import chdir, environ
    from os.path import exists, join
    from subprocess import Popen, PIPE
    try:
        install_dir = environ['DESIMODEL']
    except KeyError:
        if desimodel is not None:
            install_dir = desimodel
        else:
            install_dir = default_install_dir()
    if exists(join(install_dir, 'data')):
        raise ValueError("{0} already exists!".format(join(install_dir,
                                                           'data')))
    assert_svn_exists()

    chdir(install_dir)
    command = svn_export(version)
    # print(' '.join(command))
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    status = proc.returncode
    if status != 0:
        raise RuntimeError(err.rstrip())


def main():
    """Entry point for the :command:`install_desimodel_data` script.

    Returns
    -------
    :class:`int`
        An integer suitable for passing to :func:`sys.exit`.
    """
    from sys import argv
    from argparse import ArgumentParser
    desc = """Install desimodel data.

This script will attempt to download and install the desimodel data/ directory.
The script will attempt to attempt to install the data in the following
locations, in order of preference:

1. :envvar:`DESIMODEL`, that is, the directory specified by the
   environment variable.
2. The value set with the -d option on the command line.
3. A directory relative to the file containing this script.  This directory
   is currently {0}.

If the data directory already exists, this script will not do anything.
""".format(default_install_dir())
    parser = ArgumentParser(description=desc, prog=argv[0])
    parser.add_argument('-d', '--desimodel', action='store', dest='desimodel',
                        metavar='DESIMODEL',
                        help=('Place the data/ directory in this directory. ' +
                              'In other words, the environment variable ' +
                              'DESIMODEL should be set to this directory.'))
    parser.add_argument('-D', '--desimodel-version', action='store',
                        dest='desimodel_version', metavar='VERSION',
                        help='Explicitly set the version to download.')
    options = parser.parse_args()
    try:
        install(options.desimodel, options.desimodel_version)
    except (ValueError, RuntimeError) as e:
        print(e.message)
        return 1
    return 0
