from distutils import setup
from distutils import Extension
from Cython.Distutils import build_ext

ext = Extension("c_msdl", sources=["c_msdl.pyx"])

setup(ext_modules=[ext],
      cmdclass={'build_ext': build_ext})