import os, sys, subprocess
import os.path as osp
# import shutil


# # clean previous build
# for root, dirs, files in os.walk(".", topdown=False):
#     for name in files:
#         if (name.startswith("myext") and not(name.endswith(".pyx") or name.endswith(".pxd"))):
#             os.remove(os.path.join(root, name))
#     for name in dirs:
#         if (name == "build"):
#             shutil.rmtree(name)

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

MPI_DIR = os.environ.get('HDF5_DIR')

ext_modules=[
    cythonize(Extension("ibsimu",
        sources = ["ibsimu.pyx","src/*.cpp"],
        # library_dirs = ['.'],
        libraries = ['ibsimu'],
        depends = ['numpy>1.7'])
        )]

# Extension("kernel",
#     sources=["src/ibsimu.cpp"],
#     libraries=["gnumath"],
#     language="c++",),
#     extra_compile_args=["-fopenmp", "-O3"],
#     extra_link_args=["-DSOME_DEFINE_OPT",
#                            "-L./some/extra/dependency/dir/"]),

setup(
    name = "pybsimu",
    version = "0.01",
    author = "Pooria Heidary",
    author_email = "",
    description = ("python wrappers around ibsimu"),
    # long_description=open('README.md').read(),
    classifiers=[
                 "Development Status :: 3 - Beta",
    #             "Topic :: Utilities",
                 ],
    packages = ['pybsimu',
    # 'pybsimu.vis', 'pybsimu.formats'
    # 'pybsimu.csg', 'pybsimu.bio',
    ],
    # py_modules = modules,
    # ext_package='tracpy', 
    scripts = [],
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

# opencascade path and dependency