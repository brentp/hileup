from distutils.core import setup
from distutils.extension import Extension
import Cython
from Cython.Build import cythonize
import pysam
import numpy as np
import sys

from Cython.Compiler.Options import get_directive_defaults
version = "0.0.1"

#directive_defaults = get_directive_defaults()
#directive_defaults['linetrace'] = True
#directive_defaults['binding'] = True

setup(name="chileup",
      cmdclass={'build_ext': Cython.Build.build_ext},
      version=version,
      ext_modules = cythonize([Extension("chileup",
                               sources=["hile.c", "chileup.pyx"],
                               depends=["hile.h", "khash.h"],
                               language="c",
                               #libraries=["z", "hts"],
                               compiler_directives={'language_level', sys.version_info[0]},
                               #define_macros=[('CYTHON_TRACE', '1')],
                               extra_compile_args=['-std=c99'],
                               extra_link_args=pysam.get_libraries(),
                               include_dirs=["."] + pysam.get_include() + [np.get_include()],
                               )])
)
