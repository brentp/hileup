from distutils.core import setup
from distutils.extension import Extension
import Cython
from Cython.Build import cythonize

setup(
      cmdclass={'build_ext': Cython.Build.build_ext},
      ext_modules = cythonize([Extension("chileup",
                               sources=["hileup.c", "chileup.pyx"],
                               depends=["hileup.h", "khash.h"],
                               language='c',
                               libraries=["z", "hts"],
                               include_dirs=["../pysam/pysam",".",
                                   "/usr/include", "/usr/local/include"],
                               )])
)
