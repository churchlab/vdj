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

"""params.py

Define directory and file names that must be manually modified
to point to certain resources.
"""

import os
import warnings

warnings.simplefilter('always')

def parse_config_file(path=os.path.expanduser('~/.vdjconfig')):
    config_data = {}
    ip = open(path,'r')
    for line in ip:
        if line.startswith('#') or line.strip() == '': continue
        data = map(lambda s: s.strip(),line.split('\t'))
        config_data[data[0]] = data[1]
    ip.close()
    return config_data

config_data = parse_config_file()

# locate some directories
try:
    vdj_dir = config_data['vdj_dir']
except KeyError:
    print "Could not successfully set vdj_dir.  Does ~/.vdjconfig exist?"
    raise

try:
    imgt_dir = config_data['imgt_dir']
except KeyError:
    warning.warn("Could not find imgt_dir in .vdjconfig. May cause problems loading refseq.")

# define some other directories and variables
data_dir = 'data'
processed_dir = 'processed'

# define organism of refseq data
organism = 'human'