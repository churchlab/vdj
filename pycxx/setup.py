import os, sys
from distutils.core import setup, Extension

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

setup (name = "HelloWorld",
       packages = ['CXX'],
       package_dir = {'CXX': '.'},
       ext_modules = [
         Extension('hello',
                   sources = ['hello.cxx',
                         os.path.join(support_dir,'cxxsupport.cxx'),
                         os.path.join(support_dir,'cxx_extensions.cxx'),
                         os.path.join(support_dir,'IndirectPythonInterface.cxx'),
                         os.path.join(support_dir,'cxxextensions.c')
                         ],
            )
       ]
)
