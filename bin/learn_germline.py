import subprocess

import vdj
import vdj.analysis

# 1. Split V and J and only put in unique sequences.  Sort by lineage-weighted abundance

v_counts = pyutils.nesteddict()
j_counts = pyutils.nesteddict()
for chain in vdj.parse_imgt():
    v_feature_list = [chain.__getattribute__('V-REGION').qualifiers['gene'][0],chain.v_seq]
    v_counts.nested_add(v_feature_list)
    
    j_feature_list = [chain.__getattribute__('J-REGION').qualifiers['gene'][0],chain.j_seq]
    j_counts.nested_add(j_feature_list)

for tup in v_counts.walk():
    (keylist,val) = (tup[:-1],tup[-1])
    v_counts.nested_assign(keylist,len(val))

for tup in j_counts.walk():
    (keylist,val) = (tup[:-1],tup[-1])
    j_counts.nested_assign(keylist,len(val))

for key in v_counts:
    outhandle = open()
    for (v_seq,count) in sorted(v_counts[key].iteritems(),key=lambda t: t[1],reverse=True):
        outhandle.write(">%s\n%s\n" % (,v_seq))
    outhandle.close()


# 2. UCLUST each file and pull up to top two clusters

# start parallelization here
cmd = 'usearch --cluster %s --usersort --id %f --uc %s'
for inputfile in glob.glob("*.fasta"):
    p = subprocess.Popen(cmd % (inputfile,0.97,outputfile),shell=True)

for clusterfile in glob.glob('*.uc'):
    clusters = {}
    for line in clusterfile:
        if line[0] == 'S':
            data = line.split('\t')
            clusters.setdefault(data[-2],[]).append(data[-2])
        elif line[0] == 'H':
            data = line.split('\t')
            clusters.setdefault(data[-1],[]).append(data[-2])
    for cluster in sorted(clusters.itervalues(),key=lambda x: len(x),reverse=True)[0:2]:
        for seq in cluster:
            print >>outhandle, ">%s\n%s" % (,seq)

# end parallelization here...move to its own script and dispatch via LSF

# 3. For each cluster, generate PSSM using MUSCLE
# parallelize for LSF

# perform MUSCLE aln
for inputfilename in glob.glob('*'):
    cline = MuscleCommandline(input=infile,out=muscFile,maxiters=2,quiet='quiet',diags='diags')
    os.system(str(cline))

# generate PSSM
for inputfilename in glob.glob('*'):
    # Use Numpy to make it fast/easy

# output is 5 (a,g,t,c,-) x len of alignment, weighted by the number of unique lineages

# 4. Derive variant tree



# 5. Analyze tree to find representative seqs

# 6. Build CDR3 consenses with left/right alignment

# 7. Align new germlines to old germline and transfer annotation