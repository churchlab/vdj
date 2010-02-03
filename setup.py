# setup.py
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

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

hashcoreext = Extension(
                        "hashcore",
                        ["hashcore.c"],
                        include_dirs = get_numpy_include_dirs()
                        )

setup(  name = "vdj",
        version = "0.1",
        ext_modules = [alignmentcoreext,clusteringcoreext,hashcoreext]
    )
