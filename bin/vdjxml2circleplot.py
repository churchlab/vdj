#! /usr/bin/env python

import optparse

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm

import vdj
import vdj.analysis
import scale

option_parser = optparse.OptionParser()
option_parser.add_option('-q','--quantify',choices=['clone','junction','read'])
(options,args) = option_parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outname = args[1]
else:
    raise ValueError, "Must give an input and output argument"

features = ['v','j']
(uniq_feature_values,countdict) = vdj.analysis.vdjxml2countdict(inhandle,features,count=options.quantify)

_jet_data =   {'red':   ((0., 0, 0),
                         (0.35, 0, 0),
                         (0.66, 1, 1),
                         (0.89,1, 1),
                         (1, 0.5, 0.5)),
               'green': ((0., 0, 0),
                         (0.125,0, 0),
                         (0.375,1, 1),
                         (0.64,1, 1),
                         (0.91,0,0),
                         (1, 0, 0)),
               'blue':  ((0., 0.5, 0.5),
                         (0.11, 1, 1),
                         (0.34, 1, 1),
                         (0.65,0, 0),
                         (1, 0, 0))}


# compute the properties for the circle collection
min_count = min([min(cd.itervalues()) for cd in countdict.itervalues()])
max_count = max([max(cd.itervalues()) for cd in countdict.itervalues()])
size_scale = scale.linear(min_count,max_count).range(3,100)
color_scale = lambda c: mpl.cm.jet(scale.log(min_count,max_count).range(0,0.85)(c))

xy = []
s = []
c = []
for (i,v_gene) in enumerate(uniq_feature_values['v']):
    for (j,j_gene) in enumerate(uniq_feature_values['j']):
        try:
            count = countdict[v_gene][j_gene]
        except KeyError: # count == 0
            continue
        
        xy.append( (i,j) )
        s.append(size_scale(count))
        c.append(color_scale(count))

# make the plots
fig = plt.figure(figsize=(24,12))
ax = fig.add_subplot(111,axisbg='k')
ax.set_aspect('equal')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
ax.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
# ax.spines['bottom'].set_position(('outward',5))
# ax.spines['left'].set_position(('outward',5))
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('left')
ax.set_xlim([-2,len(uniq_feature_values['v'])+1])
ax.set_ylim([-2,len(uniq_feature_values['j'])+1])

circles = mpl.collections.CircleCollection(s,offsets=xy,transOffset=ax.transData,facecolors=c,linewidths=0,clip_on=False)

ax.add_collection(circles)
ax.autoscale_view()

# ax.set_xlabel('V')
# ax.set_ylabel('J')
# leg = ax.legend(loc=0,numpoints=1,prop=mpl.font_manager.FontProperties(size='small'))
# leg.get_frame().set_visible(False)
# fig.show()
fig.savefig(outname,bbox_inches='tight')
