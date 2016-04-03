from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ext = Extension("debriujn",
                sources = ["debruijn.pyx", "pyencoder.pyx", "pyeulertour.pyx", "component.pyx",
                           "gpu_debruijn.cu", "encoder.cu", "eulertour.cu", "component.cu"],
                language = "c"
                )


setup(
    name = "debruijn graph",
    ext_modules = cythonize([ext])
)