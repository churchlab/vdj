import tempfile
import subprocess
import time

import vdj

# ===================
# = LSF Dispatching =
# ===================

def split_into_parts(rep,outputname,packetsize):
	parts = []
	nfiles = int(math.ceil(float(len(rep)) / packetsize))
	for i in xrange(nfiles):
		# split into multiple parts; append .partNNN to filename
		currfilename = outputname + '.part%04i' % i
		parts.append(currfilename)
		vdj.writeVDJ(rep[i*packetsize:(i+1)*packetsize],currfilename)
	return parts

def split_into_good_VJCDR3s(rep,outputname,verbose=True):
	repgood = rep.get_chains_fullVJCDR3()
	parts = []
	vjcombo = []
	partnum = 0
	for vseg in refseq.IGHV[1:]:
		for jseg in refseq.IGHJ[1:]:
			currfilename = outputname + '.part%04i' % partnum
			currrep = repgood.get_chains_AND([vseg,jseg])
			if len(currrep) == 0:
				continue
			vdj.writeVDJ(currrep,currfilename,verbose=verbose)
			parts.append(currfilename)
			vjcombo.append(vseg.replace('/','_')+'_'+jseg.replace('/','_'))
			partnum += 1
	return (parts,vjcombo)

def load_parts(parts,verbose=True):
	rep = vdj.Repertoire()
	for part in parts:
		rep_part = vdj.fastreadVDJ(part,mode='Repertoire',verbose=verbose)
		rep += rep_part
	return rep

def generate_script(operation):
	scriptname = tempfile.mktemp('.py','vdj_operation','./')
	op = open(scriptname,'w')
	print >>op, "import vdj"
	print >>op, "import sys"
	print >>op, "rep = vdj.fastreadVDJ(sys.argv[1],mode='Repertoire')"
	print >>op, "rep=vdj."+operation+"(rep)"
	print >>op, "vdj.writeVDJ(rep,sys.argv[1])"
	op.close()
	return scriptname

def submit_to_LSF(queue,LSFopfile,script,parts):
	processes = []	# list of PIDs that are dispatched through LSF
	for chunk in parts:
		proc = subprocess.Popen( ['bsub','-q'+queue,'-o'+LSFopfile,'python',script,chunk], stdout=subprocess.PIPE )
		proc.wait()
		processes.append(proc.stdout.read().split('<')[1].split('>')[0])
	return processes

def waitforLSFjobs(PIDs,interval=30):
	finished = False
	while not finished:
		time.sleep(interval)
		p = subprocess.Popen('bjobs',stdout=subprocess.PIPE)
		p.wait()
		status = p.stdout.read().split('\n')
		if status[0].split()[0] != 'JOBID':
			finished = False
			continue
		runningprocesses = [line.split()[0] for line in status if line.split() != [] and line.split()[0] != 'JOBID']
		finished = True
		for pid in PIDs:
			if pid in runningprocesses:
				finished = False
	return