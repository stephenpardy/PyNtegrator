from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefiles = ['pyorbits.pyx']
ext_modules = [Extension("pyorbits",
                         sourcefiles
                         )]

setup(
    name = 'Converter',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

