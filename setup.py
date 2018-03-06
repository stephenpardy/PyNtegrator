from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

sourcefiles_pyorbit = ['src/pyorbits.pyx', 'src/halo.c', 'src/plummer.c']
sourcefiles_probe = ['src/orbit_probe.pyx']

ext_modules = [Extension("pyorbits",
                         sources=sourcefiles_pyorbit,
                         library_dirs=['/usr/local/lib', '/usr/lib'],
                         libraries=['m', 'gsl', 'gslcblas'],
                         extra_compile_args=['-std=c99', '-Wall', '-Wextra'],
                         include_dirs=[numpy.get_include()]
                         ),
                Extension("orbit_probe",
                          sources=sourcefiles_probe,
                          )]

setup(
    name = 'Converter',
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(ext_modules),
    include_dirs=[numpy.get_include()]
)
