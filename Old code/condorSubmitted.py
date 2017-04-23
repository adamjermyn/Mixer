import subprocess
import numpy as np

data = '/data/vault/asj42/'

def submit(omega, tw, w, n2):
	postfix = '_'+str(omega)+'_'+str(tw)+'_'+str(w)
	fname = 'job' + postfix
	fi = open(data + fname,'w+')
	fi.write('Universe = vanilla\n')
	fi.write('Executable = /home/asj42/Documents/Turbulence/main\n')
	fi.write('arguments= '+str(omega)+' '+str(tw)+' '+str(w)+' '+str(n2)+'\n')
	fi.write('Genenv = true\n')
	fi.write('Request_memory = 1000\n')
	fi.write('Output = '+data+'job'+postfix+'.out\n')
	fi.write('Error = '+data+'job'+postfix+'.error\n')
	fi.write('Log = '+data+'job'+postfix+'.log\n')
	fi.write('Queue\n')
	fi.close()
	subprocess.call(['condor_submit',data + fname])

omega = 10**np.linspace(-3,3,num=8,endpoint=True)
tw = np.linspace(0,2*np.pi,num=8,endpoint=True)
w = 10**np.linspace(-3,3,num=8,endpoint=True)

n2 = -1

for om in omega:
	for tww in tw:
		for ww in w:
			submit(om,tww,ww,n2)
