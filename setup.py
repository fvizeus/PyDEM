try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize

import numpy

extensions = [Extension('cDEM', ['cDEM.pyx'], include_dirs = [numpy.get_include()])]

setup(ext_modules=cythonize(extensions))