# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
============
desiutil.git
============

This package contains code for interacting with DESI git products.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# The line above will help with 2to3 support.


def last_tag(owner, repo):
    """Scan GitHub tags and return the most recent tag.

    Parameters
    ----------
    owner : :class:`str`
        The owner or group in GitHub
    repo : :class:`str`
        Name of the product.

    Returns
    -------
    :class:`str`
        The most recent tag found on GitHub.
    """
    from os.path import basename
    import requests
    api_url = 'https://api.github.com/repos/{0}/{1}/git/refs/tags/'
    r = requests.get(api_url.format(owner, repo))
    data = r.json()
    try:
        return basename(data[-1]['ref'])
    except KeyError:
        return '0.0.0'


def version(git='git'):
    """Use ``git describe`` to generate a version string.

    Parameters
    ----------
    git : :class:`str`, optional
        Path to the git executable, if not in :envvar:`PATH`.

    Returns
    -------
    :class:`str`
        A :pep:`386`-compatible version string.

    Notes
    -----
    The version string should be compatible with :pep:`386` and
    :pep:`440`.
    """
    from subprocess import Popen, PIPE
    myversion = '0.0.1.dev0'
    try:
        p = Popen([git, "describe", "--tags", "--dirty", "--always"],
                  universal_newlines=True, stdout=PIPE, stderr=PIPE)
    except OSError:
        return myversion
    except EnvironmentError:
        return myversion
    out, err = p.communicate()
    if p.returncode != 0:
        return myversion
    ver = out.rstrip().split('-')[0]+'.dev'
    try:
        p = Popen([git, "rev-list", "--count", "HEAD"],
                  universal_newlines=True, stdout=PIPE, stderr=PIPE)
    except OSError:
        return myversion
    except EnvironmentError:
        return myversion
    out, err = p.communicate()
    if p.returncode != 0:
        return myversion
    ver += out.rstrip()
    return ver
