from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefiles = ['pyorbits.pyx']
ext_modules = [Extension("pyorbits",
                         sources=sourcefiles,
                         library_dirs = ['/usr/local/lib'],
                         libraries = ['m', 'gslcblas', 'gsl']
                         )]

setup(
    name = 'Converter',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

