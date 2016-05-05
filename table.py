import numpy as np
from subprocess import Popen, PIPE

def parseOut(out):
	lines = out.split('\n')
	lines = [l.split(' ') for l in lines][:-1]
	lines = [[float(l) for l in line if l!=''] for line in lines]
	return np.array(lines)

def run(s,p,om,n,lvl):
	pp = Popen(['./main','0','0','0',str(om),'0','0',str(s),str(p),str(n),'0',str(lvl)], stdin=PIPE, stdout=PIPE, stderr=PIPE)
	output, err = pp.communicate()
	rc = pp.returncode
	d = None
	try:
		d = parseOut(output)
	except Exception as e:
		print 'Error: Switching to debug version!'
		pp = Popen(['./debug','0','0','0',str(om),'0','0',str(s),str(p),str(n),'0',str(lvl)], stdin=PIPE, stdout=PIPE, stderr=PIPE)
		output, err = pp.communicate()
		rc = pp.returncode
		try:
			d = parseOut(output)
		except Exception as e:
			print output,rc
			print '-----'
			print i,j,k,l
			print '-----'
			print  ' '.join(['./main','0','0','0',str(om),'0','0',str(s),str(p),str(n),'0',str(lvl)])
			print '-----'
	return d

def lvlFinder(s,p,om,n,stlvl=3,rtol=1e-3,atol=1e-3):
	level = stlvl
	while True:
		d0 = run(s,p,om,n,level)
		if d0 is None:
			return None
		d1 = run(s,p,om,n,level+1)
		err = d0 - d1
		err = np.reshape(err,(-1,))
		err = np.abs(err)
		err = (err < rtol*np.reshape(np.abs(d1),(-1,))+atol)
		if np.sum(err) == len(err):
			return d1
		level += 1
		print level
		if level > 15:
			print 'Could not satisfy error criterion.'
			return d1


