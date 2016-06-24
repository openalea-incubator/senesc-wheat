# -*- coding: latin-1 -*-
"""
Notes:

- use setup.py develop when tracking in-development code
- when removing modules or data files from the project, run setup.py clean --all and delete any obsolete .pyc or .pyo.

"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import ez_setup
ez_setup.use_setuptools()

import sys
from setuptools import setup, find_packages

import senescwheat

if sys.version_info < (2, 7):
    print('ERROR: Senesc-Wheat requires at least Python 2.7 to run.')
    sys.exit(1)

if sys.version_info >= (3, 0):
    print('WARNING: Senesc-Wheat has not been tested with Python 3.')

setup(
    name = "Senesc-Wheat",
    version=senescwheat.__version__,
    packages = find_packages(),

    include_package_data = True,

    # metadata for upload to PyPI
    author = "C.Chambon, R.Barillot",
    author_email = "camille.chambon@grignon.inra.fr, romain.barillot@grignon.inra.fr",
    description = "Model of senescence adapted to wheat ",
    long_description = "Mod�le de senescence adaptee au bl�",
    license = "", # TODO
    keywords = "", # TODO
    url = "https://sourcesup.renater.fr/projects/senesc-wheat/",
    download_url = "", # TODO
)
