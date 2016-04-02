from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ext = Extension("debriujn",
                sources = ["debruijn.pyx", "encoder.pyx", "eulertour.pyx", "component.pyx", 
                           "gpu_debruijn.cu", "gpu_encoder.cu", "gpu_eulertour.cu", "gpu_component.cu"],
                language = "c"
                )


setup(
    name = "debruijn graph",
    ext_modules = cythonize([ext])
)