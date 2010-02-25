#!/bin/bash
# pycxx_ghetto_setup.sh
# 
# This script installs PyCXX for python 2.x.
# 
# To install PyCXX in an inelegant way, all that needs to happen is to copy some
# relevant source files to the proper locations. Then build your extensions
# using distutils and make sure to point to the correct source files for PyCXX.
# 
# The script pycxx_ghetto_setup.sh should copy all the right files to the proper
# locations. It's pretty easy to undo what it does by looking at the script.
# 
# Note: it uses `python` as the executable to obtain some basic information
# about the install If you use some other name for the executable, change it
# (eg., 'python25').

# obtain directory of Python installation
PREFIX=`python -c "import sys; print sys.prefix"`
PY_VERSION_SHORT=`python -c "import sys; print str(sys.version_info[0])+'.'+str(sys.version_info[1])"`

# eliminate references to Python2 directory in the code
mv CXX/Python2/Extensions.hxx CXX/Python2/Extensions.hxx.original
sed "s/\/Python2//g" CXX/Python2/Extensions.hxx.original > CXX/Python2/Extensions.hxx

# create necessary directories
mkdir -p $PREFIX/share/python$PY_VERSION_SHORT/CXX
mkdir -p $PREFIX/include/python$PY_VERSION_SHORT/CXX
mkdir -p $PREFIX/lib/python$PY_VERSION_SHORT/site-packages/CXX

# copy files to correct locations
cp Src/Python2/*    $PREFIX/share/python$PY_VERSION_SHORT/CXX
cp CXX/Python2/*    $PREFIX/include/python$PY_VERSION_SHORT/CXX
cp CXX/WrapPython.h $PREFIX/include/python$PY_VERSION_SHORT/CXX
cp CXX/Version.hxx  $PREFIX/include/python$PY_VERSION_SHORT/CXX
cp Lib/__init__.py  $PREFIX/lib/python$PY_VERSION_SHORT/site-packages/CXX