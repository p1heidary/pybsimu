import os, sys, subprocess
import os.path as osp
from glob import glob
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

if glob('*.o') and 'clean' in sys.argv:
    [os.remove(f) for f in glob('*.o')]

MPI_DIR      = os.environ.get('MPI_DIR')
IBSIMU_DIR   = os.environ.get('IBSIMU_DIR')
# source_files = glob(os.path.join('src', '*.cpp'))
# source_files.insert(0, "kernel.pyx")

ext_modules = cythonize(Extension(
        "kernel",
        sources = ["kernel.pyx"],
        language="c++",
        include_dirs = ['src', '/usr/include/cairo',
            '/usr/include/gtk-2.0', '/usr/include/glib-2.0',
            '/usr/lib64/glib-2.0/include', # /usr/bin/glib-config
            '/usr/include/pango-1.0', # /usr/lib64/pkgconfig/pango.pc
            '/usr/lib64/gtk-2.0/include/', # /usr/lib64/pkgconfig/gdk-2.0.pc
            '/usr/include/gdk-pixbuf-2.0', # /usr/lib64/pkgconfig/gdk-pixbuf-2.0.pc
            '/usr/include/atk-1.0', # /usr/lib64/pkgconfig/atk.pc
            ], # need config script
        libraries=[ IBSIMU_DIR + "/lib64/libibsimu-1.0.6.so"],
        ))

# Extension("kernel",
#     sources=["src/ibsimu.cpp"],
#     libraries=["gnumath"],
#     language="c++",),
#     depends = ['numpy>1.7 ']
#     library_dirs = ['.'],
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
    # ext_package='trcpy', 
    scripts = [],
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

# opencascade path and dependency
# python setup.py build_ext --inplace