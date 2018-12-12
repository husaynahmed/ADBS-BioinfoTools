##############################################
##############################################
#############    ViV ver 1.0    ##############
############# Venn in VCFs 1.0  ##############
##############################################
##############################################

#######################################################################
### ViV is a part of the ADBS-ToolKit(adbs-tk) and is a tool 
### to find the number of common genetic variants across multiple 
### samples. This tool can simultaneously process thousands 
### of VCF files obtained from whole genome/exome sequencing.
### Please read the README file before using this program.
#######################################################################

# packages import
import pandas as pd
import numpy as np
import argparse



###
##### Argument for dendrogen
###

parser = argparse.ArgumentParser(description='ViV version 1.0', epilog='Please give absolute(full) path to all the files. Refer the README at https://github.com/husaynahmed/adbs-tk for more information.')

parser.add_argument('-I', '--inputfileinfo',
                     help='Path to the text file containing sample IDs and list of VCF file paths. One entry per line.',
                     required='True')

parser.add_argument('-O', '--outdir',
                     help='Path to the output directory where the results will be written.',
                     required='True')

args = parser.parse_args()
inputFI = args.inputfileinfo
outdir = args.outdir
#####################################################################

###
##### Reading input files
###
input = pd.read_table(inputFI, header=None, prefix="C",comment='#', usecols=(0,1))
print "||| Reading input files information |||"
sam_names = input["C0"].tolist()
vcf_files = input["C1"].tolist()

len1 = len(sam_names)
len2 = len(vcf_files)
print "You entered the following samples", sam_names

##### Opening the vcf files and storing in variables

def frNameCr(nd):
	"Create the names of dataframes"
	frName = []
	for i in range(0,nd):
		frName.append("df_"+str(i))
	return frName

print "||| Reading VCF files |||"
j = 0
frP = frNameCr(len1)
for i in vcf_files:
	df = pd.read_table(i, header=None, prefix="C",comment='#', usecols=(0,1,3,4))
	frP[j] = df	
	j=j+1
#####################################################################

###
##### Identifying common variants across samples and creating the Venn matrix
###

print "||| Performing ViV |||"
j = 0
x= []
venn_mat = np.array(x, dtype=float)

for i in sam_names:
	k = 0
	dfA = frP[j]
	for h in range(0,len(frP)):
		dfB = frP[k]
		mer = pd.merge(dfA, dfB, how='inner')
		print len(mer),"\t",
		venn_mat = np.append(venn_mat, len(mer))
		k=k+1
	j=j+1


venn_matA = np.reshape(venn_mat, (len1,len1))
#####################################################################


###
##### Writing the results in a matrix form to 
###
print "||| Writing results to output directory |||"
df_matC = pd.DataFrame(venn_matA, columns=sam_names, index=sam_names)
out_file = outdir + "/Venn_in_VCFs_results.tsv"
df_matC.to_csv(out_file, mode='a', sep='\t')

print "||| ViV - COMPLETED !!! |||"
#####################################################################



#######################################################################################
###  Developed: Husayn Ahmed P 
###  https://github.com/husaynahmed
###  At ADBS, National Centre for Biological Sciences, NCBS, Bengaluru
###
#######################################################################################

