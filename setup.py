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
    'maligner',
    sources = ['maligner.cxx'],
    )


setup(  name = "vdj",
        version = "1.2",
        packages = ['CXX'],
        package_dir = {'CXX': '.'},
        include_dirs= [r'.','/usr/include/python2.6','/usr/include/python2.6/CXX'],
        library_dirs= [r'.'],
        libraryes=['stdc++','m'],
        ext_modules = [maligner]
    )
