#### Python
# calculate how significantly each annotation improves the fit

import os
from scipy import stats
import re
import argparse

cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Calculate pvalue of annotations improvement to PAINTOR model fit')
parser.add_argument('-l', '--locus', nargs='?', default='locus', help='locus name to be included in the output name')
parser.add_argument('-w', '--wrkdir', nargs='?', default=cwd, help='path to directory containing the PAINTOR output files to test')
parser.add_argument('-o', '--outdir', nargs='?', default=cwd, help='path to directory for output')
parser.add_argument('-b', '--baseName', nargs='?', default="BF.Base", help='name of the bayes factor file for the base model (no annotations)')
parser.add_argument('-f', '--outputname', nargs='?', default="", help='filename for output file containing annotation p-values')

args = parser.parse_args()

locus = args.locus
wd = args.wrkdir
od = args.outdir
baseFileName = args.baseName
outputname = args.outputname


# wd = '/media/BurchardRaid01/LabShare/Home/pgoddard/wrkdir_lungfxn_gwas_sage/results/12_paintor_finemap_admixmap/locus1'
# od = '/media/BurchardRaid01/LabShare/Home/pgoddard/wrkdir_lungfxn_gwas_sage/results/12_paintor_finemap_admixmap'

def ratioTest(logBFbase, logBFmodel):
	lr = -2*(float(logBFbase) - float(logBFmodel))
	pval = 1-stats.chi2.cdf(lr, 1)
	return pval

def getPval(outFileName, baseFile, prefix = "BF."):
	'''
	baseFile = path to your PAINTOR base model Bayes Factor results (default is BF.Base)
	modelDir = path to the directory containing the locus-phenotype-specific Bayes Factor results for all the tested annotations
	prefix = prefix for the bayes factor files (default = "BF.")
	'''

	# get base model
	with open(baseFile, 'r') as Base:
		logBFbase = Base.readline().split()[0]
		logBFbase = float(logBFbase)

	with open(outFileName, "w") as outFile:

		# write header
		outFile.write("annotation\tlog(BF)\tpval\n")

		# get list of annotation results
		for file in os.listdir('.'):
			fileName = os.path.basename(file)

			if fileName == baseFile:
				continue

			# calculate pvalue for each annotation output
			elif fileName.startswith(prefix):
				annot = re.sub('BF.locus\d+.', '', fileName)
				with open(fileName, 'r') as inFile:

					# get the model and base bayes factors
					logBFmodel = inFile.readline().split()[0]
					logBFmodel = float(logBFmodel)

					# perform the ratio test
					pval = ratioTest(logBFbase, logBFmodel)

					# save results
					data = [annot, str(logBFmodel), str(pval)]
					outFile.write('\t'.join(data) + '\n')


def main():
	print('Processing ' + locus)
	os.chdir(wd)
	#outName = od + '/' + locus + '.annotations.pvals'
	outName = od + '/' + outputname
	print 'results will be written to ' + outName
	getPval(outFileName = outName, baseFile = baseFileName)

if __name__== "__main__":
  main()
