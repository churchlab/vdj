# setup.py
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

import os
#os.environ['CXX'] = 'g++'
os.environ['CPP'] = 'g++'
#os.environ['LDSHARED'] = 'g++'

alignmentcoreext = Extension(
                        "alignmentcore",
                        ["alignmentcore.c"],
                        include_dirs = get_numpy_include_dirs()
                        )

alignmentcore2ext = Extension(
                        "alignmentcore2",
                        ["alignmentcore2.cpp"],
                        include_dirs = get_numpy_include_dirs()
                        )

clusteringcoreext = Extension(
                        "clusteringcore",
                        ["clusteringcore.c"],
                        include_dirs = get_numpy_include_dirs()
                        )

setup(  name = "vdj",
        version = "0.1",
        ext_modules = [alignmentcoreext,clusteringcoreext,alignmentcore2ext]
    )
