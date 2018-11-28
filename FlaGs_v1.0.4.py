#Author: Chayan Kumar Saha, Gemma C. Atkinson
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import random
import argparse
import ftplib
import os, sys, os.path, math
import gzip
import getopt
from collections import OrderedDict
from tkinter import *
master = Tk()
import subprocess
Entrez.email = "gemma.atkinson@gmail.com"

usage= ''' Description:  Identify flanking genes and cluster them based on similarity and visualize the structure; Requirement= Python3, BioPython; tkinter ; Optional Requirement= ETE3'''


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-a", "--assemblyList", help="Protein Accession with assembly Identifier eg. GCF_000001765.3 in a text file separated by newline")
parser.add_argument("-p", "--proteinList", help="Protein Accession eg. XP_ or WP_047256880.1 in a text file separated by newline")
parser.add_argument("-r", "--redundant",action="store_true", help="Search all GCFs for each query")
parser.add_argument("-e", "--ethreshold", help="E value threshold. Default = 1e-10")
parser.add_argument("-n", "--number", help="Number of Jackhmmer iterations. Default = 3")
parser.add_argument("-s", "--strand", help="Number of strand for looking up or downstream. Default = 4")
parser.add_argument("-t", "--tree", action="store_true", help="If you want to see flanking genes along with phylogenetic tree, requires ETE3 installation. By default it will not produce.")
parser.add_argument("-ts", "--tshape", help="Size of triangle shapes that represent flanking genes, this option only works when -t is used. Default = 12 ")
parser.add_argument("-tf", "--tfontsize", help="Size of font inside triangles that represent flanking genes, this option only works when -t is used. Default = 4 ")
parser.add_argument("-to", "--tree_order", action="store_true", help="Generate Output with Tree, and then use the tree order to generate other view. ")
parser.add_argument("-o", "--out_prefix", required= True, help="Any Keyword to define your output eg. MyQuery")
parser.add_argument("-k", "--keep", action="store_true", help="If you want to keep the intermediate files eg. gff3 use [-k]. By default it will remove.")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.3')
parser.add_argument("-vb", "--verbose", action="store_true", help="Use this option to see the work progress for each query as stdout. ")
args = parser.parse_args()
parser.parse_args()

if args.assemblyList:
	if args.redundant:
		print('"-r" option is only works with "-p", please try again with proper command')
		sys.exit()

if args.tree:
	if args.tshape:
		if int(args.tshape)>0:
			size=int(args.tshape)
		else:
			print("Kindly input the size of triangles, recommended 12. Not applicable for 0 and negative values")
	else:
		size=12
else:
	if args.tshape:
		print("Kindly make sure that you are using -t to make this -ts argument working")
		sys.exit()

if args.tree:
	if args.tfontsize:
		if int(args.tfontsize)>0:
			fsize=str(args.tfontsize)
		else:
			print("Kindly input the font Size required inside triangles, recommended 4. Not applicable for 0 and negative values")
	else:
		fsize=str(4)
else:
	if args.tfontsize:
		print("Kindly make sure that you are using -t to make this -tf argument working")
		sys.exit()

if args.ethreshold:
	evthresh=args.ethreshold
else:
	evthresh="1e-10"

if args.number:
	iters=args.number
else:
	iters="3"

if args.strand:
	if int(args.strand)>0:
		s= str(int(args.strand)+1)
	else:
		print('Please insert positive values, starting from 1')
		sys.exit()
else:
	s= "5"

queryList=[] #Preparation of format input as a list.


if args.proteinList and not args.assemblyList:
	with open (args.proteinList, 'r') as pList:
		for line in pList:
			if line[0]!='#':
				Line=line.rstrip().split('\t')
				if len(Line)>1:
					print('Check Input file, Incorrect Format.')
					sys.exit()
				else:
					queryList.append(Line)
elif args.assemblyList and not args.proteinList :
	with open (args.assemblyList, 'r') as apList:
		for line in apList:
			if line[0]!='#':
				Line=line.rstrip().split('\t')
				if len(Line)<2:
					print('Check Input file, Incorrect Format.')
					sys.exit()
				else:
					newFormat=Line[1]+'\t'+Line[0]
					queryList.append(newFormat.split('\t'))
else:
	print('Incorrect Input!')
	sys.exit()

def accession_from_xp(accession_nr):
	"""
	:param accession_nr: NCBI protein accession
	:return: Bioproject number of all species for that protein which is used to grab Assembly number
	"""
	#Entrez.email = "gemma.atkinson@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
	# blocked, warning is sent to this email.
	try:
		handle = Entrez.efetch(db="protein", id=accession_nr, rettype="gbwithparts", retmode="text")
	except Exception as e:
		print(str(e) + ", error in entrez-fetch protein accession, %s, not found in database. \n"
		"Continuing with the next protein in the list. \n" % accession_nr)
		return False
	record = SeqIO.read(handle, "genbank")
	bioproj=record.dbxrefs
	handle.close()
	bio=[]
	for item in bioproj:
		if item.split(':')[0]=='BioProject':
			bio.append(item.split(':')[1])
	return set(bio)

def accession_from_wp(accession_nr):

	"""
	:param accession_nr: NCBI protein accession
	:return: Set of assembly number of all species for particular protein

	"""
	#Entrez.email = "gemma.atkinson@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
	# blocked, warning is sent to this email.
	try:
		handle = Entrez.efetch(db="protein", id=accession_nr, rettype="ipg", retmode="xml")
	except Exception as e:
		print(str(e), ", error in entrez-fetch protein accession, {}, not found in database. \n" "Continuing with the next protein in the list. \nError in function: {}".format(accession_nr, accession_from_wp.__name__))
		return False
	record = list(Entrez.parse(handle))
	handle.close()
	if len(record)>0:
		if "ProteinList" in record[0]:
			for keys in (record[0]["ProteinList"]):
				assembly=[]
				if "CDSList" in keys:
					for item in keys["CDSList"]:
						if "assembly" in item.attributes :
							assembly.append(item.attributes["assembly"])
					return set(assembly)
				else:
					return("NAI")
		else:
			return("NAI")
	else:
		return("NAI")

def seq_from_wp(accession_nr):
	"""
	:param accession_nr: NCBI protein accession
	:return: Protein Sequence
	"""
	if accession_nr[-1]!='*':
		#Entrez.email = "gemma.atkinson@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
		# blocked, warning is sent to this email.
		try:
			handle = Entrez.efetch(db="protein", id=accession_nr, rettype="gbwithparts", retmode="text")
		except Exception as e:
			print(str(e), ", error in entrez-fetch protein accession, {}, not found in database. \n" "Continuing with the next protein in the list. \nError in function: {}".format(accession_nr, accession_from_wp.__name__))
			return False
		record = SeqIO.read(handle, "genbank")
		handle.close()
		return record.description.split('[')[0]+'\t'+record.seq
	else:
		return accession_nr[:-1]+'\t'+'--'

def des_check(item):
	if item:
		return item
	else:
		return 'notFound'


def normalize_strand(item1, item2):  #Strand change
	if item1=='+':
		return item2
	else:
		if item2=='+':
			return '-'
		else:
			return '+'

def up(item):
	if item=='+':
		return 'Upstream '
	else:
		return 'Downstream '

def down(item):
	if item=='+':
		return 'Downstream '
	else:
		return 'Upstream '

def ups(item):
	if item=='+':
		return '-'
	else:
		return '+'

def downs(item):
	if item=='+':
		return '+'
	else:
		return '-'


ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
ftp.cwd("/genomes/refseq") # move to refseq directory

filenames = ftp.nlst() # get file/directory names within the directory
if 'assembly_summary_refseq.txt' in filenames:
	ftp.retrbinary('RETR ' + 'assembly_summary_refseq.txt', open('input.txt', 'wb').write) # get the assembly summary from refseq

assemblyName={}
bioDict={} #bioproject as keys and assemble number (eg.GCF_000001765.1) as value
accnr_list_dict={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
with open('input.txt', 'r') as fileIn:
	for line in fileIn:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			accnr_list_dict[Line[0]]= Line[7]+'\t'+Line[19]
			bioDict[Line[1]]=Line[0]
			assemblyName[Line[0]]=Line[0]

ftp_gen = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
ftp_gen.cwd("/genomes/genbank") # move to refseq directory

filenames = ftp_gen.nlst() # get file/directory names within the directory
if 'assembly_summary_genbank.txt' in filenames:
	ftp_gen.retrbinary('RETR ' + 'assembly_summary_genbank.txt', open('input_genbank.txt', 'wb').write) # get the assembly summary from refseq


assemblyName_GCA={}
bioDict_gen={}
accnr_list_dict_gen={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
with open('input_genbank.txt', 'r') as fileIn:
	for line in fileIn:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			if len(Line)>19:
				if Line[18]=='identical':
					if Line[17] in accnr_list_dict:
						bioDict_gen[Line[1]]=Line[0]
						accnr_list_dict_gen[Line[0]]= accnr_list_dict[Line[17]]
						assemblyName_GCA[Line[0]]=Line[17]

bioDict.update(bioDict_gen)
accnr_list_dict.update(accnr_list_dict_gen)
assemblyName.update(assemblyName_GCA)

print ('\n'+ '>> Database Downloaded. Work in progress ...'+ '\n')




q=0
ne=0
queryDict={} #protein Id as query and a set of assembly number as value [either All or Species of interest]
with open (args.out_prefix+'_NameError.txt', 'w') as fbad:
	for query in queryList:
		q+=1
		if len(query)<2:
			if query[0][0]!='X' and query[0][1:3]=='P_':
				queryDict[query[0]+'#'+str(q)]=accession_from_wp(query[0])
			elif query[0][:2]=='XP':
				assemList=[]
				for bioprojs in accession_from_xp(query[0]):
					if bioprojs in bioDict:
						assemList.append(bioDict[bioprojs])
				queryDict[query[0]+'#'+str(q)]=set(assemList)
			else:
				ne+=1
				print(query[0], file= fbad)

		else:
			if query[0][0]!='X' and query[0][1:3]=='P_':
				asset=set()
				if len(accession_from_wp(query[0]))>0:
					for elements in accession_from_wp(query[0]):
						if query[1]==elements:
							asset.add(query[1])
				queryDict[query[0]+'#'+str(q)]=asset
			elif query[0][:2]=='XP':
				asset=set()
				for bioprojs in accession_from_xp(query[0]):
					if bioDict[bioprojs]==query[1]:
						asset.add(query[1])
				queryDict[query[0]+'#'+str(q)]=asset
			else:
				ne+=1
				print(query[0], file= fbad)

#print(queryDict)
nai=0
NqueryDict={}
with open (args.out_prefix+'_Insufficient_Info_In_DB.txt', 'w') as fNai:
	for query in queryDict:
		if queryDict[query]:
		#if len(queryDict[query])!=0:
			if args.redundant:
				NqueryDict[query]=queryDict[query]
			else:
				NqueryDict[query]=random.sample(queryDict[query],1)
		else:
			print(query, file=fNai)
			nai+=1
newQ=0
for query in NqueryDict:
	newQ+=1

print('>> Input file assessment report: ')
print('\t'+'Discarded protein ids with improper accession prefix : '+str(ne)+'. See "'+args.out_prefix+'_NameError.txt'+'" file for details.')
print('\t'+'Discarded protein ids lacking proper information in RefSeq DB : '+str(nai)+'. See "'+args.out_prefix+'_Insufficient_Info_In_DB.txt'+'" file for details.')
print('\t'+'Remaining queries: '+str(newQ))

FoundDict={} #Accession that found in Refseq
FlankFoundDict={} #Accession that have flanking genes
accFlankDict={} #{'WP_092250023.1#1': {0: 'WP_092250023.1+', 1: 'WP_092250020.1+', 2: 'WP_092250017.1-', -1: 'tRNA*+', -2: 'WP_092250026.1-'}}
positionDict={} #Accession as keys:Start and end position as value
speciesDict={} #SpeciesName stored here
queryStrand={} #Strand Information for each query
LengthDict={} #Length of each query

seqDict={}
desDict={}

acc_CGF_Dict={}

#with open(args.out_prefix+'_flankgene.tsv', 'w') as ftab:
count=0
for query in NqueryDict:
	count+=1
	if args.verbose:
		print('\n'+'> '+str(count)+' in process out of '+str(newQ)+' ... '+'\n')
		print('Query Name =', query.split('#')[0], '\n')
	for item in NqueryDict[query]:
		a=0
		LineList=[]
		geneProt={} # 'gene2504': 'WP_092248795.1', 'gene1943': 'tRNA'
		geneChrom={} #'gene1708': 'NZ_MJLP01000030.1'
		#item='GCF_'+items[4:]
		if item in accnr_list_dict:
			ftpLine=accnr_list_dict[item].split('\t')[1]
			ftpsplitDir = ftpLine.split('/')[3:]
			ftp_path = '/'.join(map(str,ftpsplitDir))
			ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
			ftp.cwd('/'+ftp_path)
			files = ftp.nlst()
			for elements in files:
				if 'gff.gz' in elements: #Check if GFF.gz is there
					#time.sleep(0.1)
					try:
						ftp.retrbinary('RETR ' + elements, open(assemblyName[item]+'.gff.gz', 'wb').write)
						with gzip.open(assemblyName[item]+'.gff.gz', 'rb') as gffIn: #Download and read gff.gz
							for line in gffIn:
								if line.decode('utf-8')[0]!='#':
									Line=line.decode('utf-8').rstrip().split('\t')
									if Line[2]=='CDS':
										if Line[8].split(';')[3][:5]=='Name=': #eliminates pseudo gene as they don't have 'Name='
											geneProt[Line[8].split(';')[1].split('=')[1]]=Line[8].split(';')[3].split('=')[1]
											geneChrom[Line[8].split(';')[1].split('=')[1]]=Line[0]
									if Line[2]=='gene':
										a+=1
										newGene=str(a)+'\t'+Line[8].split(';')[0][3:]+'\t'+ Line[3]+'\t'+Line[4]+'\t'+ Line[6]+ '\t'+ Line[0]
										LineList.append(newGene.split('\t')) #1       gene3006        10266   10342   -       NZ_FPCC01000034.1
										for genDes in Line[8].split(';'):
											if 'gene_biotype=' in genDes:
												if Line[8].split(';')[0][3:] not in geneProt:
													geneProt[Line[8].split(';')[0][3:]]=genDes.split('=')[1]+'_'+query.split('#')[1]+'*'
							geneList=[] ##List of gene names coding same protein (accession), we are taking one from them
							for genes in geneProt:
								if query.split('#')[0]==geneProt[genes]:
									geneList.append(genes)
							rangeList=[]
							for line in LineList:
								if geneChrom[geneList[0]]==line[5]:
									rangeList.append(int(LineList[LineList.index(line)][0]))
							for genes in geneProt:
								if genes==geneList[0]:
									if query.split('#')[0]==geneProt[genes]:
										for line in LineList:
											if genes==line[1]:
												FoundDict[query]='Yes'
												speciesDict[query]=accnr_list_dict[item].split('\t')[0]
												queryStrand[query]= LineList[LineList.index(line)][4]
												#print(LineList[LineList.index(line)])
												positionDict[query]= ("\t".join(map(str,LineList[LineList.index(line)][2:-2])))
												LengthDict[query]= int(LineList[LineList.index(line)][3])-int(LineList[LineList.index(line)][2])+1
												udsDict={}
												dsDict={}
												udsDict[0]= query+'+' # O strand
												for x in range(1,int(s)):
													if int(LineList[LineList.index(line)-x][0]) in rangeList:
														acc_CGF_Dict[query]= LineList[LineList.index(line)-x][-1] +'\t'+ item +'\t'+ ftpLine
														seqDes=str(seq_from_wp(geneProt[LineList[LineList.index(line)-x][1]]))
														seqDict[geneProt[LineList[LineList.index(line)-x][1]]]=seqDes.split('\t')[1]
														desDict[geneProt[LineList[LineList.index(line)-x][1]]]=seqDes.split('\t')[0]
														positionDict[geneProt[LineList[LineList.index(line)-x][1]]+'#'+query.split('#')[1]]= ("\t".join(map(str,LineList[LineList.index(line)-x][2:-2])))
														LengthDict[geneProt[LineList[LineList.index(line)-x][1]]+'#'+query.split('#')[1]]= int(LineList[LineList.index(line)-x][3])-int(LineList[LineList.index(line)-x][2])+1
														udsDict[int(ups(LineList[LineList.index(line)][4])[0]+str(x))]= geneProt[LineList[LineList.index(line)-x][1]]+'#'+query.split('#')[1]+\
															normalize_strand(LineList[LineList.index(line)][4],LineList[LineList.index(line)-x][4])
												for y in range(1,int(s)):
													if int(LineList[LineList.index(line)+y][0]) in rangeList:
														acc_CGF_Dict[query]= LineList[LineList.index(line)+y][-1] +'\t'+ item +'\t'+ ftpLine
														seqDes=str(seq_from_wp(geneProt[LineList[LineList.index(line)+y][1]]))
														seqDict[geneProt[LineList[LineList.index(line)+y][1]]]=seqDes.split('\t')[1]
														desDict[geneProt[LineList[LineList.index(line)+y][1]]]=seqDes.split('\t')[0]
														positionDict[geneProt[LineList[LineList.index(line)+y][1]]+'#'+query.split('#')[1]]= ("\t".join(map(str,LineList[LineList.index(line)+y][2:-2])))
														LengthDict[geneProt[LineList[LineList.index(line)+y][1]]+'#'+query.split('#')[1]]= int(LineList[LineList.index(line)+y][3])-int(LineList[LineList.index(line)+y][2])+1
														dsDict[int(downs(LineList[LineList.index(line)][4])[0]+str(y))]= geneProt[LineList[LineList.index(line)+y][1]]+'#'+query.split('#')[1]+\
															normalize_strand(LineList[LineList.index(line)][4],LineList[LineList.index(line)+y][4])
												udsDict.update(dsDict)
												accFlankDict[query]=udsDict
												if query in accFlankDict:
													if len(accFlankDict[query])>0:
														FlankFoundDict[query]='Yes'
														if args.verbose:
															print('\t', query.split('#')[0], 'Report: Flanking Genes Found', '\n')
													else:
														FlankFoundDict[query]='No'
														if args.verbose:
															print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')
												else:
													FlankFoundDict[query]='No'
													if args.verbose:
														print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')
							gff_gz='./'+assemblyName[item]+'.gff.gz'
							if args.keep:
								pass
							else:
								if os.path.isfile(gff_gz):
									os.remove(gff_gz)
								else:    ## Show an error ##
									print("Error: %s file not found" % gff_gz)
					except:
						FlankFoundDict[query]='No'
						FoundDict[query]='No: ProteinID was not found in Genome Assembly'
						pass
						if args.verbose:
							print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')

						#print(query.split('#')[0], item, ftp_path, sep='\t', file=error)
		else:
			FoundDict[query]='Yes : Assembly ID did not match with RefSeq '
			pass
			if args.verbose:
				print('\t', 'Corresponding Assembly ID', item, 'does not exist in NCBI Refseq', '\n')

#accFlankDict
'''
with open (args.out_prefix+'.log', 'w') as logout:
#print(accFlankDict)
	print('queryDict=',queryDict, file=logout)
	print('NqueryDict=',NqueryDict, file=logout)
	print('assemblyName=', assemblyName, file=logout)
	print('FoundDict=',FoundDict, file=logout)#Accession that found in Refseq
	print('FlankFoundDict=',FlankFoundDict, file=logout)#Accession that have flanking genes
	print('accFlankDict=',accFlankDict, file=logout) #{'WP_092250023.1#1': {0: 'WP_092250023.1+', 1: 'WP_092250020.1+', 2: 'WP_092250017.1-', -1: 'tRNA*+', -2: 'WP_092250026.1-'}}
	print('positionDict=',positionDict, file=logout)#Accession as keys:Start and end position as value
	print('speciesDict=',speciesDict, file=logout)#SpeciesName stored here
	print('queryStrand=',queryStrand, file=logout)#Strand Information for each query
	print('LengthDict=',LengthDict, file=logout) #Length of each query
	print('seqDict=',seqDict, file=logout)
	print('desDict=',desDict, file=logout)
	print('acc_CGF_Dict=',acc_CGF_Dict, file=logout)
'''

'''
{'WP_000004567.1#1': {0: 'WP_000004567.1#1+', -1: 'WP_000931496.1#1+', -2: 'WP_003310314.1#1+', -3: 'WP_033795253.1#1-', -4: 'WP_000153877.1#1+', 1: 'WP_000342211.1#1-', 2: 'WP_000708631.1#1-', 3: 'WP_000291520.1#1-', 4: 'WP_000572582.1#1-'}}
'''
allFlankGeneList=[]
for keys in accFlankDict:
	for item in accFlankDict[keys]:
		allFlankGeneList.append(accFlankDict[keys][item].split('#')[0])
#print(allFlankGeneList)

myfile="./input.txt" # for deletion of the downloaded file from ftp
if args.keep:
	pass
else:
	if os.path.isfile(myfile):
		os.remove(myfile)
	else:    ## Show an error ##
		print("Error: %s file not found" % myfile)

myfile="./input_genbank.txt" # for deletion of the downloaded file from ftp
if args.keep:
	pass
else:
	if os.path.isfile(myfile):
		os.remove(myfile)
	else:    ## Show an error ##
		print("Error: %s file not found" % myfile)

if args.keep:
	pass
else:
	subprocess.Popen("rm G*F*.gz", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)



flankF=0
with open (args.out_prefix+'_flankgene_Report.log', 'w') as errOut:
	serial=0
	print('#Serial','Query','Assembly_Found', 'FlankingGene_Found', sep='\t', file=errOut)
	for queries in NqueryDict:
		serial+=1
		if queries in FoundDict:
			if queries in FlankFoundDict:
				flankF+=1
				print(str(serial), queries.split('#')[0], FoundDict[queries], FlankFoundDict[queries], sep='\t', file=errOut)
			else:
				print(str(serial), queries.split('#')[0], FoundDict[queries], 'No', sep='\t', file=errOut)
		else:
			print(str(serial), queries.split('#')[0], 'No', 'No', sep='\t', file=errOut)

print('\n')

print('>> '+str(flankF)+' accessions ran successfully out of remaining '+str(serial)+'. See "'+args.out_prefix+'_flankgene_Report.log'+'" file for details.'+'\n'+'\n')

querySeqDict={}#{Accession_query:seq}
if args.tree:
	with open(args.out_prefix+'_tree.fasta', 'w') as treeOut:
		for queries in NqueryDict:
			if queries in FlankFoundDict:
				if FlankFoundDict[queries]=='Yes':
					handle = Entrez.efetch(db="protein", id=queries.split('#')[0], rettype="gbwithparts", retmode="text")
					record = SeqIO.read(handle, "genbank")
					if record.id==queries.split('#')[0]:
						record.id= queries+'_'+speciesDict[queries].replace(' ','_').replace(':','_').replace('[','_').replace(']','_')
						record.description= ''
						print(record.format("fasta"), file=treeOut)
						querySeqDict[str(record.id)]=str(record.seq)
					handle.close()

#print(querySeqDict)
###Going to change
#seqDict={}
#desDict={}
#with open(args.out_prefix+'_flankgene.tsv',"r") as tsvIn :
#	for line in tsvIn:
#		if line!='#':
#			line=line.rstrip().split('\t')
#			if len(line)==11:
#				desDict[line[4]]=line[9]
#				seqDict[line[4]]=line[10]
###
#print(seqDict)
#print(desDict)
if len(seqDict)!=len(desDict):
	if len(seqDict)>len(desDict):
		for seqids in sorted(seqDict):
			if seqids not in desDict:
				desDict[seqids]=des_check(str(seq_from_wp(seqids).split('\t')[0]))
	else:
		for seqids in sorted(desDict):
			if seqids not in seqDict:
				seqDict[seqids]=str(seq_from_wp(seqids).split('\t')[1])
else:
	if args.verbose:
		print ('Description collected for Flanking Genes!')


b=0
with open (args.out_prefix+'_flankgene.fasta'+'_cluster_out', 'w') as fastaNew:
	for seqids in sorted(seqDict):
		if seqDict[seqids]!='--':
			b+=1
			print('>'+seqids+'|'+desDict[seqids]+'\n'+seqDict[seqids], file=fastaNew)

if args.verbose:
	print ('Total Flanking genes found = '+ str(b))



infilename=args.out_prefix+'_flankgene.fasta'+'_cluster_out'

directory = args.out_prefix+'_flankgene.fasta'+'_cluster_out_individuals'
if not os.path.exists(directory):
    os.makedirs(directory)

infile=open(args.out_prefix+'_flankgene.fasta'+'_cluster_out',"r").read()
al=infilename+"_"+iters+"_"+evthresh+"_jackhits.tsv"
outacclists=open(al,"w")


i=1
for seqids in sorted(seqDict):
	if seqDict[seqids]!='--':
		i_f=directory+"/"+str(i)+".txt"
		indivfile=open(i_f,"w")
		indivfile.write(">"+seqids+'\n'+seqDict[seqids])
		indivfile.close()
		command="jackhmmer -N %s --incE %s --incdomE %s --tblout %s/tblout%s.txt %s  %s>%s/out%s.txt" %(iters, evthresh, evthresh, directory, str(i), i_f, infilename, directory, str(i))
		#print(command)
		os.system(command)
		tbl=open(directory+"/tblout"+str(i)+".txt").read()
		part=tbl.split("----------\n")[1].split("\n#")[0]
		lines=part.splitlines()
		acclist=[]
		for line in lines:
			inc=line.split()[17]
			#acc=line.split("|")[-3].split(" ")[0]
			acc=line.split('|')[0]
			if inc=="1":
				acclist.append(acc)
		outacclists.write(str(i)+"\t"+str(acclist)+"\n")
		i=i+1

outacclists.close()

raw=open(al).read().strip()

d={}
index=0
for line in raw.split("\n"):
	if line.split("\t")[1]!='[]':
		index+=1
		actxt=line.split("\t")[1].replace(",","").replace("[","").replace("]","").replace("'","")
		actlist=actxt.split(" ")
		d[index]=(actlist)


i=1
while i<len(d)+1:	#use i and j to iterate through the combinations
    list1=d[i]
    j=i+1
    while j<len(d)+1:
    	list2=d[j]
		#print "i", i, list1, " vs ",
		#print "j", j, list2
    	if set(list1) & (set(list2)): # if there is an intersection
			#print "yes there is", list1, list2
    		union=list(set(list2) | set(list1))
    		d[j]=union
    		d[i]=[] #...and empty the redundant list
    	j=j+1
    i=i+1


#print("==")
#d : {1: ['WP_000153877.1'], 2: [], 3: ['WP_000291520.1'], 4: ['WP_000342211.1'], 5: [],... 30: ['WP_055032751.1', 'WP_001246052.1']}

trueAccessionCount={}
for keys in d:
	numbers=[]
	for item in d[keys]:
		numbers.append(allFlankGeneList.count(item))
	trueAccessionCount[(';'.join(map(str,d[keys])))]=sum(numbers)

#trueAccessionCount: 'WP_001229260.1;WP_001229255.1': 4
odtrue=OrderedDict(sorted(trueAccessionCount.items(), key= lambda item:item[1],reverse=True))
#print(odtrue)


familyNumber=0
with open(infilename+"_"+iters+"_"+evthresh+"_clusters.tsv","w") as clusOut:
	for k, v in odtrue.items():
		#print (k.split(';'), odtrue[k], v)
		if len(k.split(';'))>0 and v>0:
			familyNumber+=1
			print(str(familyNumber),str(odtrue[k]),k, sep='\t', file=clusOut)
#			k=k+1


outfile_des=open(infilename+"_"+iters+"_"+evthresh+"_outdesc.txt","w")
inf=open(infilename+"_"+iters+"_"+evthresh+"_clusters.tsv","r")

acclists=inf.read().splitlines()
for line in acclists:
	acclist=line.split("\t")[2].split(";")
	familyAssignedValue=line.split("\t")[0]
	if int(line.split("\t")[1])>1:
		for acc in acclist:
            #print (acc, "\t", desDict[acc])
			outfile_des.write(familyAssignedValue+'('+str(allFlankGeneList.count(acc))+')'+"\t"+acc+"\t"+desDict[acc]+"\n")
        #print ("\n\n")
		outfile_des.write ("\n\n")


import random
from random import randint
import colorsys

def random_color(h=None):
	"""Generates a random color in RGB format."""
	if not h:
		c = random.random()
	d = 0.5
	e = 0.5
	return _hls2hex(c, d, e)

def _hls2hex(c, d, e):
	return '#%02x%02x%02x' %tuple(map(lambda f: int(f*255),colorsys.hls_to_rgb(c, d, e)))

familyDict={} # Accession:Assigned family Number from Jackhammer
with open(args.out_prefix+'_flankgene.fasta_cluster_out_'+iters+'_'+evthresh+'_clusters.tsv', 'r') as clusterIn:
	for line in clusterIn:
		if line[0]!='#':
			line=line.rstrip().split('\t')
			if int(line[1])>1:
				for item in (line[2].split(';')):
					familyDict[item]=int(line[0])
			else:
				familyDict[line[2]]=0


familynum=[]
for acc in familyDict:
	familynum.append(familyDict[acc])

center=int(max(familynum))+1
noProt=int(max(familynum))+2

for ids in LengthDict:
	if ids.split('#')[0][-1]=='*':
		familyDict[ids.split('#')[0]]=noProt
	if ids.split('#')[0] not in familyDict:
		familyDict[ids.split('#')[0]]=center

color={}
color[center]='#000000'
color[noProt]='#f2f2f2'

colorDict={} #Assigned family Number from Jackhammer : colorcode
for families in set(familynum):
	if families == 0:
		colorDict[families]=str('#ffffff')
	else:
		if random_color()!='#ffffff' or random_color()!='#000000' or random_color()!='#f2f2f2' :
			colorDict[families]=random_color()


def outliner (item):
	if item =='#ffffff':
		return '#bebebe'
	elif item =='#f2f2f2':
		return '#008000'
	else:
		return item

colorDict.update(color)

maxs=(int(s)-1)
mins=maxs-(maxs*2)


if not args.tree:
	nPos=[]
	pPos=[]
	with open(args.out_prefix+'_operon.tsv', 'w') as opOut:
		for queries in accFlankDict:
			for items in sorted(accFlankDict[queries]):
				if queryStrand[queries]=='+':
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'_'+speciesDict[queries].replace(' ','_').replace(':','_').replace('[','_').replace(']','_')
					qStrand=queryStrand[queries]
					nStrand=accFlankDict[queries][items][-1]
					family=familyDict[accFlankDict[queries][items][:-1].split('#')[0]]
					startPos=int(positionDict[accFlankDict[queries][0][:-1]].split('\t')[0])-1
					start=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[0])
					end=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[1])
					if queries in acc_CGF_Dict:
						info=acc_CGF_Dict[queries]
					else:
						info='not_found'+'\t'+'not_found'+'\t'+'not_found'
					print(species, lengths, qStrand, nStrand, family, start-startPos, end-startPos, start, end, ids, info, sep='\t', file=opOut)
					nP=start-startPos
					pP=end-startPos
					nPos.append(nP)
					pPos.append(pP)

				else:
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'_'+speciesDict[queries].replace(' ','_').replace(':','_').replace('[','_').replace(']','_')
					qStrand=queryStrand[queries]
					nStrand=accFlankDict[queries][items][-1]
					family=familyDict[accFlankDict[queries][items][:-1].split('#')[0]]
					startPos=startPos=int(positionDict[accFlankDict[queries][0][:-1]].split('\t')[1])+1
					start=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[1])
					end=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[0])
					if queries in acc_CGF_Dict:
						info=acc_CGF_Dict[queries]
					else:
						info='not_found'+'\t'+'not_found'+'\t'+'not_found'
					print(species, lengths, qStrand, nStrand, family, startPos-start, startPos-end, end, start, ids, info, sep='\t', file=opOut)
					nP=startPos-start
					pP=startPos-end
					nPos.append(nP)
					pPos.append(pP)

			print('\n\n', file=opOut)

	windowMost=round(((max(ptPos)+abs(min(ntPos))+1)*4)/100)
	widthM=(windowMost*3)+500
	heightM=int(newQ)*20
	#aheightM=heightM*1.3

	canvas = Canvas(master, width=widthM,height=heightM,background='white', scrollregion=(0,0,widthM*2.5,heightM*2.5))
	hbar=Scrollbar(master,orient=HORIZONTAL)
	hbar.pack(side=BOTTOM,fill=X)
	hbar.config(command=canvas.xview)
	vbar=Scrollbar(master,orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
	#canvas.config(width=1500,height=1000)
	canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)

	def operonFamily(item):
		if item==0:
			return ' '
		elif item==center:
			return ' '
		elif item==noProt:
			return ' '
		else:
			return item


	eg1=open(args.out_prefix+'_operon.tsv','r').read()

	egs=eg1.split("\n\n\n\n")
	line_pos_y=0
	for eg in egs:
		if eg!='':
			coln=0
			entries=eg.splitlines()
			ndoms=len(entries)
			ptnstats=entries[0].split("\t")
			org=ptnstats[0].replace("_"," ")
			textspace=widthM/2
			line_pos_y=line_pos_y+16
			half_dom_height=5
			text = canvas.create_text(textspace/2,line_pos_y, text=org, fill="#404040", font=("Arial", "12"))
			for entry in entries:
				items=entry.split("\t")
				aln_start=round(int(items[5])*4/100)
				aln_end=round(int(items[6])*4/100)
				strandType=items[3]
				dom1_name=int(items[4])
				dom1_len=(aln_end-aln_start)
				oL80=round(dom1_len*80/100)
				dom1_start=aln_start+textspace
				dom1_end=dom1_len+dom1_start
				if strandType=='+':
					rect = canvas.create_polygon(dom1_start, line_pos_y+half_dom_height, dom1_start, line_pos_y-half_dom_height,dom1_start+oL80, line_pos_y-half_dom_height, dom1_end, line_pos_y, dom1_start+oL80, line_pos_y+half_dom_height,fill=colorDict[dom1_name], outline=outliner(colorDict[dom1_name]))
				else:
					rect = canvas.create_polygon(dom1_end-oL80, line_pos_y+half_dom_height, dom1_start, line_pos_y, dom1_end-oL80, line_pos_y-half_dom_height,dom1_end, line_pos_y-half_dom_height, dom1_end, line_pos_y+half_dom_height, fill=colorDict[dom1_name], outline=outliner(colorDict[dom1_name]))
				textd1 = canvas.create_text(dom1_start+(dom1_len/2),line_pos_y, text=operonFamily(dom1_name), font=("Arial", "7"))
				coln=coln+1

	retval = canvas.postscript(file=args.out_prefix+"_flankgenes.ps", height=heightM, width=widthM, colormode="color")
	mainloop()


if args.tree:###Tree Command###
	tree_file= args.out_prefix+'_tree.fasta'
	tree_command=tree_command="ete3 build -a %s -o %s --clearall -w mafft_default-trimal01-none-fasttree_full" %(tree_file, tree_file[:-6])
	#print(tree_command)
	os.system(tree_command)
	from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node

	def normalize_strandView(item):  #Strand view change
		if item=='+':
			return '>'
		else:
			return '<'

	def familyView(item):  #Strand view change
		if item==0:
			return ' '
		elif item==center:
			return ' '
		elif item==noProt:
			return ' '
		else:
			return str(item)

	seqMult=((maxs)*2)+1
	seq = ("XXXXXXXXXXXXX--"*seqMult)
	startDict={}
	udList=[]
	for ud in range (mins, maxs+1, 1):
		udList.append(ud)
	sList=[]
	for sa in range(1, 15*seqMult, 15):
		sList.append(sa)
	for ln in range(len(udList)):
		startDict[udList[ln]]=sList[ln]

	#WP_071976325.1#5_Bacillus_oleronius:0.266349
	#print(startDict)

	nwTree=''
	motifDict={}
	motifDict_2={}
	#treeOrderList=[]
	with open(args.out_prefix+'_tree/mafft_default-trimal01-none-fasttree_full/'+args.out_prefix+'_tree.fasta.final_tree.nw', 'r') as treeIn:
		for line in treeIn:
			nwTree=line
			for items in line.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
				item=items.split('_')[0]+'_'+items.split('_')[1]
				#treeOrderList.append(item)
				simple_motifs=[]
				simple_motifs_2=[]
				for keys in sorted(startDict):
					if keys in accFlankDict[item]:
						#print(keys,item, accFlankDict[item] )# accession: accFlankDict[item]['U_4'][:-1]
						simple_motifs_s = [startDict[keys], startDict[keys]+13, normalize_strandView(accFlankDict[item][keys][-1]), None, size, outliner(colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]]), 'rgradient:'+colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]], "arial|"+fsize+"|black|"+familyView(familyDict[accFlankDict[item][keys][:-1].split('#')[0]])]
						simple_motifs.append(simple_motifs_s)
						simple_motifs_2_s = [startDict[keys], startDict[keys]+13, normalize_strandView(accFlankDict[item][keys][-1]), None, size, outliner(colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]]),colorDict[familyDict[accFlankDict[item][keys][:-1].split('#')[0]]], "arial|"+fsize+"|black|"]
						simple_motifs_2.append(simple_motifs_2_s)
					else:
						simple_motifs_s = [startDict[keys], startDict[keys]+13, '[]', None, size, '#eeeeee', 'rgradient:'+'#ffffff', "arial|"+fsize+"|black|"]
						simple_motifs.append(simple_motifs_s)
						simple_motifs_2_s = [startDict[keys], startDict[keys]+13, '[]', None, size, '#eeeeee', '#ffffff', "arial|"+fsize+"|black|"]
						simple_motifs_2.append(simple_motifs_2_s)
				motifDict[items[:items.index(':')]]=simple_motifs
				motifDict_2[items[:items.index(':')]]=simple_motifs_2

	def get_example_tree():
		# Create a random tree and add to each leaf a random set of motifs
		# from the original set
		t= Tree(nwTree)
		for item in nwTree.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
			seqFace = SeqMotifFace(seq, motifs=motifDict[item[:item.index(':')]], seq_format="-", gap_format="blank")
			(t & item[:item.index(':')]).add_face(seqFace, 0, "aligned")
		t.ladderize()
		return t

	def get_example_tree_2():
		# Create a random tree and add to each leaf a random set of motifs
		# from the original set
		t= Tree(nwTree)
		for item in nwTree.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
			seqFace2 = SeqMotifFace(seq, motifs=motifDict_2[item[:item.index(':')]], seq_format="-", gap_format="blank")
			(t & item[:item.index(':')]).add_face(seqFace2, 0, "aligned")
		t.ladderize()
		return t

	if __name__ == '__main__':
		t = get_example_tree()
		ts = TreeStyle()
		ts.tree_width = 300
		ts.show_branch_support = True
		if args.tree_order:
			t.write(outfile=args.out_prefix+'_ladderTree.nw')
			t.render(args.out_prefix+"_flankgenes_1.svg",tree_style=ts)
		else:
			t.show(tree_style=ts)
			t.render(args.out_prefix+"_flankgenes_1.svg",tree_style=ts)

	if __name__ == '__main__':
		t = get_example_tree_2()
		ts = TreeStyle()
		ts.tree_width = 300
		ts.show_branch_support = True
		if args.tree_order:
			t.write(outfile=args.out_prefix+'_ladderTree.nw')
			t.render(args.out_prefix+"_flankgenes_2.svg", tree_style=ts)
		else:
			t.show(tree_style=ts)
			t.render(args.out_prefix+"_flankgenes_2.svg", tree_style=ts)

if args.tree and args.tree_order:
	treeOrderList=[]
	with open(args.out_prefix+'_ladderTree.nw', 'r') as laddertreeIn:
		for line in laddertreeIn:
			for items in line.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
				item=items.split('_')[0]+'_'+items.split('_')[1]
				treeOrderList.append(item)

	ntPos=[]
	ptPos=[]
	with open(args.out_prefix+'_TreeOrder_operon.tsv', 'w') as opOut:
		for queries in treeOrderList:
			for items in sorted(accFlankDict[queries]):
				if queryStrand[queries]=='+':
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'_'+speciesDict[queries].replace(' ','_').replace(':','_').replace('[','_').replace(']','_')
					qStrand=queryStrand[queries]
					nStrand=accFlankDict[queries][items][-1]
					family=familyDict[accFlankDict[queries][items][:-1].split('#')[0]]
					startPos=int(positionDict[accFlankDict[queries][0][:-1]].split('\t')[0])-1
					start=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[0])
					end=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[1])
					if queries in acc_CGF_Dict:
						info=acc_CGF_Dict[queries]
					else:
						info='not_found'+'\t'+'not_found'+'\t'+'not_found'
					#print(info)
					print(species, lengths, qStrand, nStrand, family, start-startPos, end-startPos, start, end, ids, info, sep='\t', file=opOut)
					nP=start-startPos
					pP=end-startPos
					ntPos.append(nP)
					ptPos.append(pP)

				else:
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'_'+speciesDict[queries].replace(' ','_').replace(':','_').replace('[','_').replace(']','_')
					qStrand=queryStrand[queries]
					nStrand=accFlankDict[queries][items][-1]
					family=familyDict[accFlankDict[queries][items][:-1].split('#')[0]]
					startPos=startPos=int(positionDict[accFlankDict[queries][0][:-1]].split('\t')[1])+1
					start=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[1])
					end=int(positionDict[accFlankDict[queries][items][:-1]].split('\t')[0])
					if queries in acc_CGF_Dict:
						info=acc_CGF_Dict[queries]
					else:
						info='not_found'+'\t'+'not_found'+'\t'+'not_found'
					#print(info)
					print(species, lengths, qStrand, nStrand, family, startPos-start, startPos-end, end, start, ids, info, sep='\t', file=opOut)
					nP=startPos-start
					pP=startPos-end
					ntPos.append(nP)
					ptPos.append(pP)

			print('\n\n', file=opOut)

	windowMost=round(((max(ptPos)+abs(min(ntPos))+1)*4)/100)
	widthM=(windowMost*3)+500
	heightM=int(newQ)*20
	#aheightM=heightM*1.3


	canvas = Canvas(master, width=widthM,height=heightM,background='white', scrollregion=(0,0,widthM*2.5,heightM*2.5))
	hbar=Scrollbar(master,orient=HORIZONTAL)
	hbar.pack(side=BOTTOM,fill=X)
	hbar.config(command=canvas.xview)
	vbar=Scrollbar(master,orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
	canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)

	def operonFamily(item):
		if item==0:
			return ' '
		elif item==center:
			return ' '
		elif item==noProt:
			return ' '
		else:
			return item


	eg1=open(args.out_prefix+'_TreeOrder_operon.tsv','r').read()

	egs=eg1.split("\n\n\n\n")
	line_pos_y=0
	for eg in egs:
		if eg!='':
			coln=0
			entries=eg.splitlines()
			ndoms=len(entries)
			ptnstats=entries[0].split("\t")
			org=ptnstats[0].replace("_"," ")
			textspace=widthM/2
			line_pos_y=line_pos_y+16
			half_dom_height=5
			text = canvas.create_text(textspace/2,line_pos_y, text=org, fill="#404040", font=("Arial", "12"))
			for entry in entries:
				items=entry.split("\t")
				aln_start=round(int(items[5])*4/100)
				aln_end=round(int(items[6])*4/100)
				strandType=items[3]
				dom1_name=int(items[4])
				dom1_len=(aln_end-aln_start)
				oL80=round(dom1_len*80/100)
				dom1_start=aln_start+textspace
				dom1_end=dom1_len+dom1_start
				if strandType=='+':
					rect = canvas.create_polygon(dom1_start, line_pos_y+half_dom_height, dom1_start, line_pos_y-half_dom_height,dom1_start+oL80, line_pos_y-half_dom_height, dom1_end, line_pos_y, dom1_start+oL80, line_pos_y+half_dom_height,fill=colorDict[dom1_name], outline=outliner(colorDict[dom1_name]))
					print(dom1_start, line_pos_y+half_dom_height, dom1_start, line_pos_y-half_dom_height,dom1_start+oL80, line_pos_y-half_dom_height, dom1_end, line_pos_y, dom1_start+oL80, line_pos_y+half_dom_height)
				else:
					rect = canvas.create_polygon(dom1_end-oL80, line_pos_y+half_dom_height, dom1_start, line_pos_y, dom1_end-oL80, line_pos_y-half_dom_height,dom1_end, line_pos_y-half_dom_height, dom1_end, line_pos_y+half_dom_height, fill=colorDict[dom1_name], outline=outliner(colorDict[dom1_name]))
					print(dom1_end-oL80, line_pos_y+half_dom_height, dom1_start, line_pos_y, dom1_end-oL80, line_pos_y-half_dom_height,dom1_end, line_pos_y-half_dom_height, dom1_end, line_pos_y+half_dom_height)
				textd1 = canvas.create_text(dom1_start+(dom1_len/2),line_pos_y, text=operonFamily(dom1_name), font=("Arial", "7"))
				coln=coln+1

	retval = canvas.postscript(file=args.out_prefix+"_treeOrder_flankgenes.ps", height=heightM, width=widthM, colormode="color")




print('\n'+'<<< Done >>>')
sys.exit()
