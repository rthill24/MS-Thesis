from distutils.core import setup
from Cython.Build import cythonize

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

setup(
    ext_modules = cythonize("c_msdl.pyx", annotate=True)
)