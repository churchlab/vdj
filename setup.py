# setup.py
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs
from distutils.core import setup, Extension
import os, sys

support_dir = os.path.normpath(
                   os.path.join(
			sys.prefix,
			'share',
			'python%d.%d' % (sys.version_info[0],sys.version_info[1]),
			'CXX') )

if os.name == 'posix':
	CXX_libraries = ['stdc++','m']
else:
	CXX_libraries = []

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

maligner = Extension('maligner',
                sources = ['malign.cpp', 'maligner.cxx',
                    os.path.join(support_dir,'cxxsupport.cxx'),
                    os.path.join(support_dir,'cxx_extensions.cxx'),
                    os.path.join(support_dir,'IndirectPythonInterface.cxx'),
                    os.path.join(support_dir,'cxxextensions.c')
                ],
            )

setup(  name = "vdj",
        version = "2.0",
        packages = ['CXX'],
        package_dir = {'CXX': '.'},
        include_dirs= [r'.','/usr/include/python2.6','/usr/include/python2.6/CXX'],
        library_dirs= [r'.'],
        libraryes=['stdc++','m'],
        ext_modules = [maligner, alignmentcoreext,clusteringcoreext]
    )
