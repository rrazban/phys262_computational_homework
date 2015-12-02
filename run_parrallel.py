#!/usr/local/bin/python

#run ising_MCMC.py multiple times in parralel

import os,sys,math
import numpy as np
from multiprocessing import *
import ising_MCMC
from datetime import datetime

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('grid_size', type=int, metavar='N',
                        help="Grid size; system will be NxN")
    parser.add_argument('temperature', type=float, metavar='TEMP', help="units of kT")
    parser.add_argument('--num-sweeps', type=int, metavar='NUM', default=2500)
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--out-pkl', required=False, metavar='OUT.pkl',help="pkl log file")
    parser.add_argument('--H', type=float, default=0., metavar='FIELD',
                        help="external magnetic field")
    parser.add_argument('--num-sims', type=int, metavar='NUM', default=1)
    return parser.parse_args()

def mp_run(par, nproks):
	def worker(vers):
		for ver in vers:
			dest='v'+str(ver)+'.pkl'
			ising_MCMC.main(par,dest)

	chunksize=int(math.ceil(par.num_sims/float(nproks)))
	procs=[]
	for i in range(nproks):
		p=Process(target=worker,args=(range(par.num_sims)[chunksize*i:chunksize*(i+1)],)) #the comma allows is to be read as an array and not individual elements
		procs.append(p)
		p.start()

	for p in procs:
		p.join()

if __name__=='__main__':
	nprocs=8		#number proccessors over which to parallize, make sure nprocs !> max processors of your computer
	parameters=parse_args()
	begin=str(datetime.now())
	mp_run(parameters,nprocs)
	end=str(datetime.now())

	print 'start time: '+begin
	print 'end time: '+end
