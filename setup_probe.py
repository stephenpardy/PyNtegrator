from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

sourcefiles = ['orbit_probe.pyx']
ext_modules = [Extension("orbit_probe",
                         sources=sourcefiles,
                         )]

setup(
    name = 'Converter',
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(ext_modules)
)
