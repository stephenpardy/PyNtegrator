from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

sourcefiles = ['pyorbits.pyx', 'halo.c']
ext_modules = [Extension("pyorbits",
                         sources=sourcefiles,
                         library_dirs = ['/usr/local/lib'],
                         libraries = ['m', 'gslcblas', 'gsl'],
                         extra_compile_args= ['-std=c99']
                         )]

setup(
    name = 'Converter',
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(ext_modules)
)
