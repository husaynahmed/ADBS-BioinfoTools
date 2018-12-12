##############################################
##############################################
############# dendrogen ver 1.0 ##############
##############################################
##############################################

#######################################################################
### dendrogen is a part of the ADBS-ToolKit(adbs-tk) and is built   ###
### to perform clustering of samples based on the sharing  		  ###
### of genetic variants among them. This tool can read thousands of ###
### VCFs (Variant call files) generated from whole genome or exome  ###
### sequencing and cluster them based on the allelic sharing.  	  ###
### A cluster dendrogram is generated that helps in visualizing	  ###
### cluster groups. dendrogen is developed as part of the 		  ###
### Accelerator program for Discovery in Brain disorders using 	  ###
### Stem cells (ADBS) at NCBS. 						  ###
### Please read the README file before using this program.		  ###
#######################################################################

# packages import
import pandas as pd
import numpy as np
import os
import glob
import argparse



###
##### Argument for dendrogen
###

parser = argparse.ArgumentParser(description='dendrogen version 1.0', epilog='Please give absolute(full) path to all the files. Refer the README at https://github.com/husaynahmed/adbs-tk for more information.')

parser.add_argument('-I', '--inputfileinfo',
                     help='Path to the text file containing list of VCF file paths. One VCF file per line.',
                     required='True')

parser.add_argument('-O', '--outdir',
                     help='Path to the output directory where the dendrogen results will be written.',
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

j = 0
frP = frNameCr(len1)
for i in vcf_files:
	df = pd.read_table(i, header=None, prefix="C",comment='#', usecols=(0,1,3,4))
	frP[j] = df	
	j=j+1
#####################################################################

###
##### Identifying common variants across samples
###

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
##### Calculating the distance matrix 
###

len3 = len(venn_matA)
matB = np.empty(shape=[len3, len3])

for i in range(0,len3):
	for j in range(0,len3):
		if( i == j ):
			matB[i,j] = 0
		else:
			matB[i,j] = (venn_matA[i,j] / venn_matA[j,j]) + (venn_matA[i,j] / venn_matA[i,i])

matC = np.empty(shape=[len3, len3])
for i in range(0,len3):
	for j in range(0,len3):
		if( i == j ):
			matC[i,j] = 0
		else:
			matC[i,j] = ( 1 / matB[i,j] )

df_matC = pd.DataFrame(matC, columns=sam_names, index=sam_names)
out_file = outdir + "/distanceMatrix_forPlotting.csv"
df_matC.to_csv(out_file, mode='a', sep=',')
#####################################################################


###
##### Performing clustering and plotting dendrogram (Using R)
###

os.chdir(outdir)
os.system("R -q -e \"matC = read.csv('distanceMatrix_forPlotting.csv', row.names = 1); distMat = as.dist(matC); hclust(distMat); pdf('result_dendrogram.pdf',width=11.693,height=8.2675); plot (hclust(distMat, method='complete')); dev.off();\"")
#####################################################################



#######################################################################################
###  Developed: Husayn Ahmed P 
###  https://github.com/husaynahmed
###  At ADBS, National Centre for Biological Sciences, NCBS, Bengaluru
###
###  Please acknowledge the usage of dendrogen by citing the following :
###  Suhas Ganesh,  Husayn Ahmed P,  Ravi K Nadella, Ravi P More, Manasa Sheshadri, Biju Viswanath, Mahendra Rao, Sanjeev Jain, The ADBS consortium, Odity Mukherjee. 2018. Exome sequencing in families with severe mental illness identifies novel and rare variants in genes implicated in Mendelian neuropsychiatric syndromes. Psychiatry and Clinical Neurosciences. doi:10.1111/pcn.12788
###
#######################################################################################

