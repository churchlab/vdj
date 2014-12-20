# Copyright 2014 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# setup.py
from setuptools import setup
from numpy.distutils.core import Extension
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

setup(  name = "vdj",
        version = "1.6.0",
        ext_modules = [alignmentcoreext,clusteringcoreext]
    )
