# setup.py
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

import os
os.environ['cc'] = 'g++'
os.environ['CXX'] = 'g++'
os.environ['CPP'] = 'g++'
os.environ['LDSHARED'] = 'g++'

alignmentcoreext = Extension(
                        "alignmentcore",
                        ["alignmentcore.c"],
                        include_dirs = get_numpy_include_dirs()
                        )

clusteringcoreext = Extension(
                        "clusteringcore",
                        ["clusteringcore.c"],
                        include_dirs = get_numpy_include_dirs()
                        )

maligner = Extension(
    "maligner",
    ["malign.pyx", "malign.cpp"],
    language="c++",
    include_dirs= [r'.','/usr/lib/', '/opt/local/include/'],
    library_dirs= [r'.'],
    libraries=["stdc++","stl"]
    )


setup(  name = "vdj",
        version = "1.2",
        ext_modules = [alignmentcoreext,clusteringcoreext, maligner]
    )
