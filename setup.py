from distutils.core import setup
from distutils.extension import Extension
import Cython
from Cython.Build import cythonize
import pysam

setup(
      cmdclass={'build_ext': Cython.Build.build_ext},
      ext_modules = cythonize([Extension("chileup",
                               sources=["hile.c", "chileup.pyx"],
                               depends=["hile.h", "khash.h"],
                               language='c',
                               libraries=["z", "hts"],
                               include_dirs=["."] + pysam.get_include(),
                               )])
)
