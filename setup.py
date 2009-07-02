# setup.py
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

setup( 	name = "alignmentcore",
		version = "0.1",
		ext_modules = [
			Extension(
				"alignmentcore",
				["alignmentcore.c","alignmentcorewrapper.c"],
				include_dirs = get_numpy_include_dirs()
			)
		]
	)
