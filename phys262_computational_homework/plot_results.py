#!/usr/bin/python
"""
plot_results.py

designed for pickle files of ising_MCMC.py 
"""
import sys,glob
import numpy as np
import cPickle as pkl
import matplotlib.pyplot as plt


class Plot_IsingSystem:
	def __init__(self):
		self.files=glob.glob('*.pkl')

	def get_finalME(self,ophile):
		stuff=pkl.load(ophile)
		phinalM=stuff['magnetizations'][len(stuff['magnetizations'])-1]
		phinalE=stuff['energies'][len(stuff['energies'])-1]
		return phinalM,phinalE 
		
	def get_datapoints(self):
		shape=(2,len(self.files))
		phinaldat=np.zeros(shape,dtype='float')
		for i,afile in enumerate(self.files):
			ofile=open(afile,'r')
			finalM,finalE=self.get_finalME(ofile)
			phinaldat[0][i]=finalM
			phinaldat[1][i]=finalE				
		return phinaldat
	
	def avg_datapoints(self,finaldat1):
		shape=(2,len(self.files))
		phinalavg=np.zeros(shape,dtype='float')
		for i in range(len(self.files)):
			phinalavg[0][i]=finaldat1[0][:i+1].sum()/float(i+1)	
			phinalavg[1][i]=finaldat1[1][:i+1].sum()/float(i+1)	
		return phinalavg

	def plot(self,finalavg1):
		fig=plt.figure()
		ax=fig.add_subplot(111)
		ax.plot(finalavg1[0],'-b',label='Magnetization')
		ax.legend(loc=2)
		ax2=ax.twinx()
		ax2.plot	
		ax2.plot(finalavg1[1],'-r',label='Energy')

		ax2.legend(loc=1)
		ax.set_xlabel("Monte-Carlo Time")
		plt.show()

	def driver(self):
		finaldat=self.get_datapoints()
		finalavg=self.avg_datapoints(finaldat)	
		self.plot(finalavg)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pkl-group', action="append", nargs='+')
    return parser.parse_args()

def main(args):
	print '# plot_results.py'
	plot_isingsystem=Plot_IsingSystem()
	plot_isingsystem.driver()

if __name__ == "__main__":
    main(parse_args())
