#!/usr/bin/env python3

import numpy as np
import lammps_logfile
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="calculate properties from lammps log files")
parser.add_argument("-i", "--input_file", type=str, default="log.lammps",  help="LAMMPS output file")
parser.add_argument("-p", "--properties",  nargs="+", type=str, help="Properties to calculate")
parser.add_argument("-o", "--out_file",type=str, default='out.txt', help="Out file name")
parser.add_argument("-plot", "--plot", default=False, action='store_true', help="Plot the block average")
parser.add_argument("-ps", "--plot_step", type=int, default=5, help="Plot step")
parser.add_argument("-m", "--maxBlockSize",type=int,default=0, help="Max block size")
parser.add_argument("-b", "--begin_step",type=int,default=4000, help="Consider the properties from this step")
args = parser.parse_args()

def blockAverage(datastream, isplot=False, maxBlockSize=0, figname='blockAverage.png'):
	"""This program computes the block average of a potentially correlated timeseries "x", and
	provides error bounds for the estimated mean <x>.
	As input provide a vector or timeseries "x", and the largest block size.

	Check out writeup in the following blog posts for more:
	http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty_14.html
	http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html
	"""

	Nobs         = len(datastream)           # total number of observations in datastream
	minBlockSize = 1;                        # min: 1 observation/block

	if maxBlockSize == 0:
		maxBlockSize = int(Nobs/4);        # max: 4 blocs (otherwise can't calc variance)

	NumBlocks = maxBlockSize - minBlockSize   # total number of block sizes

	blockMean = np.zeros(NumBlocks)               # mean (expect to be "nearly" constant)
	blockVar  = np.zeros(NumBlocks)               # variance associated with each blockSize
	blockCtr  = 0

				#
				#  blockSize is # observations/block
				#  run them through all the possibilities
				#

	for blockSize in range(minBlockSize, maxBlockSize):

		Nblock    = int(Nobs/blockSize)               # total number of such blocks in datastream
		obsProp   = np.zeros(Nblock)                  # container for parcelling block

		# Loop to chop datastream into blocks
		# and take average
		for i in range(1,Nblock+1):

			ibeg = (i-1) * blockSize
			iend =  ibeg + blockSize
			obsProp[i-1] = np.mean(datastream[ibeg:iend])

		blockMean[blockCtr] = np.mean(obsProp)
		blockVar[blockCtr]  = np.var(obsProp)/(Nblock - 1)
		blockCtr += 1

	v = np.arange(minBlockSize,maxBlockSize) # detailed info

	if isplot:
		s = args.plot_step
		plt.subplot(2,1,1)
		plt.plot(v[::s], np.sqrt(blockVar)[::s],'ro-',lw=2)
		plt.xlabel('block size')
		plt.ylabel('std')

		plt.subplot(2,1,2)
		plt.errorbar(v[::s], blockMean[::s], np.sqrt(blockVar)[::s])
		plt.ylabel('<x>')
		plt.xlabel('block size')

#		print "<x> = {0:f} +/- {1:f}\n".format(blockMean[-1], np.sqrt(blockVar[-1]))

		plt.tight_layout()
		plt.savefig(figname, dpi=300)
		plt.close()

	return blockMean[-1], np.sqrt(blockVar)[-1]

log = lammps_logfile.File(args.input_file)
f = args.out_file
b = args.begin_step

if args.properties:
    props = args.properties
else:
    props = ["Temp", "Press", "Lx", "Ly", "Lz", "PotEng", "KinEng", "TotEng"]

Data = dict()
new_vol = False
if 'Lx' in props and 'Ly' in props and 'Lz' in props:
	new_vol = True
	props.remove('Lx')
	props.remove('Ly')
	props.remove('Lz')
	props.append('Volume')

for i in props:
	Data[i] = dict()
	if i == 'Volume' and new_vol:
		data = log.get('Lx')[b:] * log.get('Ly')[b:] * log.get('Lz')[b:]
	else:
		data = log.get(i)[b:]
	# Data[i]['avg'] = np.mean(data)
	# Data[i]['std'] = np.std(data)
	Data[i]['block_avg'], Data[i]['block_std'] = blockAverage(data, args.plot, args.maxBlockSize, i + '.png')

with open(args.out_file, 'w') as file:
	file.write('#' + ' '.join(props) + '\n')
	file.write(' '.join([str(Data[i]['block_avg']) for i in props]) + '\n')
	file.write(' '.join([str(Data[i]['block_std']) for i in props]) + '\n')
	# file.write(' '.join([str(Data[i]['avg']) for i in props]) + '\n')
	# file.write(' '.join([str(Data[i]['std']) for i in props]) + '\n')
