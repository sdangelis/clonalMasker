#--------------------------------------------------------------------------------
# Copyright (c) 2020 Michael A. Boemo (mb915@cam.ac.uk)

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#--------------------------------------------------------------------------------


import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

VERSION="0.1.0"

#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """refine.py: Refines the boundary of CNAs called with low resolution.
To run refine.py, do:
  python refine.py --highRes /path/to/highRes.bed --lowRes /path/to/lowRes.bed -o /path/to/outputPrefix
Required arguments are:
     --highRes              input bed file with CNVs called at high resolution,
     --lowRes               input bed file with CNVs called at low resolution,
  -o,--output               output prefix for bed file.
Optional arguments are:
  -b,--buffer               extends the boundary of low resolution SCNAs when computing overlap (default: 0 bp).
Written by Michael Boemo, Department of Pathology, University of Cambridge.
Please submit bug reports to mb915@cam.ac.uk."""

	print(s)
	print('Version ',VERSION)
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.buffer = 0

	for i, argument in enumerate(args):
			
		if argument == '--highRes':
			a.highResPath = str(args[i+1])

		elif argument == '--lowRes':
			a.lowResPath = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outPrefix = str(args[i+1])

		elif argument == '-b' or argument == '--buffer':
			a.buffer = int(args[i+1])

		elif argument == '-h' or argument == '--help' or argument == '--version' or argument == '-v':
			splashHelp()

	#check that required arguments are met
	if not ( hasattr( a, 'highResPath') or hasattr( a, 'lowResPath') or hasattr( a, 'outPrefix') ):
		splashHelp() 
	return a


#-------------------------------------------------
def chr2len(chromosome):

	if chromosome == 'chr1':
		return 248956422
	elif chromosome == 'chr2':
		return 	242193529
	elif chromosome == 'chr3':
		return 	198295559	
	elif chromosome == 'chr4':
		return 	190214555	
	elif chromosome == 'chr5':
		return 	181538259	
	elif chromosome == 'chr6':
		return 	170805979	
	elif chromosome == 'chr7':
		return 	159345973
	elif chromosome == 'chr8':
		return 	145138636
	elif chromosome == 'chr9':
		return 	138394717
	elif chromosome == 'chr10':
		return 	133797422
	elif chromosome == 'chr11':
		return 	135086622
	elif chromosome == 'chr12':
		return 	133275309	
	elif chromosome == 'chr13':
		return 	114364328	
	elif chromosome == 'chr14':
		return 	107043718
	elif chromosome == 'chr15':
		return 	101991189
	elif chromosome == 'chr16':
		return 	90338345	
	elif chromosome == 'chr17':
		return 	83257441	
	elif chromosome == 'chr18':
		return 	80373285
	elif chromosome == 'chr19':
		return 	58617616
	elif chromosome == 'chr20':
		return 	64444167
	elif chromosome == 'chr21':
		return 	46709983	
	elif chromosome == 'chr22':
		return 	50818468	
	elif chromosome == 'chrX':
		return 	156040895
	else:
		print('chromosome not recognised')


#-------------------------------------------------
def chr2centromere(chromosome):

	if chromosome == 'chr1':
		return 121700000, 125100000
	elif chromosome == 'chr2':
		return 	91800000, 96000000
	elif chromosome == 'chr3':
		return 	87800000, 94000000	
	elif chromosome == 'chr4':
		return 	48200000, 51800000	
	elif chromosome == 'chr5':
		return 	46100000, 51400000	
	elif chromosome == 'chr6':
		return 	58500000, 62600000	
	elif chromosome == 'chr7':
		return 	58100000, 62100000
	elif chromosome == 'chr8':
		return 	43200000, 47200000
	elif chromosome == 'chr9':
		return 	42200000, 45500000
	elif chromosome == 'chr10':
		return 	38000000, 41600000
	elif chromosome == 'chr11':
		return 	51000000, 55800000
	elif chromosome == 'chr12':
		return 	33200000, 37800000	
	elif chromosome == 'chr13':
		return 	16500000, 18900000	
	elif chromosome == 'chr14':
		return 	16100000, 18200000
	elif chromosome == 'chr15':
		return 	17500000, 20500000
	elif chromosome == 'chr16':
		return 	35300000, 38400000	
	elif chromosome == 'chr17':
		return 	22700000, 27400000	
	elif chromosome == 'chr18':
		return 	15400000, 21500000
	elif chromosome == 'chr19':
		return 	24200000, 28100000
	elif chromosome == 'chr20':
		return 	25700000, 30400000
	elif chromosome == 'chr21':
		return 	10900000, 13000000	
	elif chromosome == 'chr22':
		return 	13700000, 17400000	
	elif chromosome == 'chrX':
		return 	58100000, 63800000
	else:
		print('chromosome not recognised')


#-------------------------------------------------
def somy2var(somy):

	splitS = somy.split('-')
	return int(splitS[0])


#-------------------------------------------------
class SCNA:
	def __init__(self, chromosome, start, end, var):
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.var = var


#-------------------------------------------------
class SingleCellSeq:
	def __init__(self):
		self.SCNAs = []
		self.name = ''


#-------------------------------------------------
def parseBed(filename):

	f = open(filename,'r')
	#print filename
	allCells = {}
	first = True
	for line in f:

		if len(line.rstrip()) == 0:
			continue

		if line[0:5] == 'track':

			if not first:
				allCells[cell.name] = cell
			first = False
			cell = SingleCellSeq()
			splitLine = line.rstrip().split('"')
			splitLine2 = splitLine[1].split(' ')
			cell.name = splitLine2[-1:][0][0:splitLine2[-1:][0].index('.')]
		else:

			splitLine = line.rstrip().split('\t')
			chromosome = splitLine[0]
			start = int(splitLine[1])
			end = int(splitLine[2])
			call = somy2var(splitLine[3])

			#print chromosome, start, end, call
			if call != 2:

				cell.SCNAs.append(SCNA(chromosome,start,end,call))
			
	allCells[cell.name] = cell
	return allCells


#-------------------------------------------------
def countRes(res):
	count = 0
	for c in range(1,23):
		for i in range(0,chr2len('chr'+str(c)),res):
			count += 1
	for i in range(0,chr2len('chrX'),res):
		count += 1
	return count


#-------------------------------------------------
def columnCmap(resolution,palette):

	row = []
	for c in range(1,23):
		for i in range(0,chr2len('chr'+str(c)),resolution):
			row.append(palette[c])
	for i in range(0,chr2len('chrX'),resolution):
		row.append(palette[c])
	return row


#-------------------------------------------------
def cellToRow(cell, resolution):

	row = []
	for c in range(1,23):
		for i in range(0,chr2len('chr'+str(c)),resolution):

			inSCNA = False
			for s in cell.SCNAs:
				if s.start < i and i < s.end and s.chromosome == 'chr'+str(c):
					row.append(s.var)
					inSCNA = True
					break
			if not inSCNA:
				row.append(2)
	return row

	
#MAIN---------------------------------------------
args = parseArguments(sys.argv[1:])

#import the cells
cells_highRes = parseBed(args.highResPath)
cells_lowRes = parseBed(args.lowResPath)
cells_refined = []

for key in cells_lowRes:

	refined_cell = SingleCellSeq()
	refined_cell.name = key
	refined_SCNAs = []

	if key in cells_highRes:
	
		for lowCNA in cells_lowRes[key].SCNAs:

			low_lb = int(lowCNA.start)
			low_ub = int(lowCNA.end)
			low_cn = int(lowCNA.var)

			refined_lb = low_lb
			refined_ub = low_ub

			max_ub = 0
			min_lb = 10000000000

			adjusted = False
			for highCNA in cells_highRes[key].SCNAs:
			
				high_lb = int(highCNA.start)
				high_ub = int(highCNA.end)
				high_cn = int(highCNA.var)

				if low_cn != high_cn or lowCNA.chromosome != highCNA.chromosome:
					continue

				overlap = max(0, min(low_ub+args.buffer,high_ub) - max(low_lb-args.buffer,high_lb))

				if overlap > 0:

					if high_lb < min_lb:
						min_lb = high_lb

					if high_ub > max_ub:
						max_ub = high_ub

					adjusted = True

			if adjusted:

				refined_lb = min_lb
				refined_ub = max_ub

			refined_SCNAs.append(SCNA(lowCNA.chromosome,refined_lb,refined_ub,low_cn))
		
		#glue together SCNAs that now overlap following refinement
		continueStitching = True

		while continueStitching:

			continueStitching = False
			stitched_SCNAs = []

			for i, s1 in enumerate(refined_SCNAs):

				stitchedThis = False
				for j,s2 in enumerate(refined_SCNAs[i+1:]):

					if s1.var != s2.var or s1.chromosome != s2.chromosome:
						continue

					overlap = max(0, min(s1.end,s2.end) - max(s1.start,s2.start))
					if overlap > 0:

						stitched_SCNAs.append(SCNA(s1.chromosome,min(s1.start,s2.start),max(s1.end,s2.end),s1.var))
						stitchedThis = True
						continueStitching = True
						del refined_SCNAs[i+1+j]
						break

				if not stitchedThis:
					stitched_SCNAs.append(s1)

			refined_SCNAs = stitched_SCNAs

		refined_cell.SCNAs = refined_SCNAs
		cells_refined.append(refined_cell)
		
	else:
		cells_refined.append(cells_lowRes[key])

#write refined bed
f_out = open(args.outPrefix + '_refined'+'.bed','w')
for cell in cells_refined:
	f_out.write('track name="refined CNVs, CNV state for '+cell.name+'.bam" description="refined CNVs, CNV state for '+cell.name+'.bam" visibility=1 itemRgb=On priority=96\n')
	for s in cell.SCNAs:
		f_out.write(str(s.chromosome)+'\t'+str(s.start)+'\t'+str(s.end)+'\t'+str(s.var)+'-somy\n')
f_out.close()

matrix = []
resolution = 2000000 #2 Mb
for cell in cells_refined:

	row = np.array(cellToRow(cell, resolution))
	matrix.append(row)

plt.figure()
colPalette = sns.color_palette("dark",n_colors=23)
colColours = columnCmap(resolution,colPalette)
p = sns.clustermap(np.array(matrix),col_colors=colColours,metric='euclidean',method='complete',row_cluster=True,col_cluster=False,cmap=sns.diverging_palette(10,240, as_cmap=True,center="light"),center=2.0,vmin=0,vmax=4)
plt.savefig(args.outPrefix + '_refined'+'.pdf')
plt.close()

matrix = []
for cell in cells_lowRes.values():

	row = np.array(cellToRow(cell, resolution))
	matrix.append(row)

plt.figure()
colPalette = sns.color_palette("dark",n_colors=23)
colColours = columnCmap(resolution,colPalette)
p = sns.clustermap(np.array(matrix),col_colors=colColours,metric='euclidean',method='complete',row_cluster=True,col_cluster=False,cmap=sns.diverging_palette(10,240, as_cmap=True,center="light"),center=2.0,vmin=0,vmax=4)
plt.savefig(args.outPrefix + '_lowRes'+'.pdf')
plt.close()

matrix = []
for cell in cells_highRes.values():

	row = np.array(cellToRow(cell, resolution))
	matrix.append(row)

plt.figure()
colPalette = sns.color_palette("dark",n_colors=23)
colColours = columnCmap(resolution,colPalette)
p = sns.clustermap(np.array(matrix),col_colors=colColours,metric='euclidean',method='complete',row_cluster=True,col_cluster=False,cmap=sns.diverging_palette(10,240, as_cmap=True,center="light"),center=2.0,vmin=0,vmax=4)
plt.savefig(args.outPrefix + '_highRes'+'.pdf')
plt.close()

