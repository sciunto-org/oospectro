from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("oospectro/third_party/_peak_finding_utils.pyx"),
)

