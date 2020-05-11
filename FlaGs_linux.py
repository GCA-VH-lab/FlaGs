__author__		= "Chayan Kumar Saha, Gemma C. Atkinson"
__copyright__	= "GNU General Public License v3.0"
__email__		= "chayan.sust7@gmail.com"

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import math, re
import argparse
import ftplib
import socket
import random
import time
from random import randint
import colorsys
import os, sys, os.path, math
import gzip
import getopt
from collections import OrderedDict
import subprocess
from tkinter import *
os.environ['DISPLAY'] = ':0' #linux
master = Tk()


usage= ''' Description:  Identify flanking genes and cluster them based on similarity and visualize the structure; Requirement= Python3, BioPython; tkinter ; Optional Requirement= ETE3. '''


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-a", "--assemblyList", help=" Protein Accession with assembly Identifier eg. GCF_000001765.3 in a text file separated by newline. ")
parser.add_argument("-p", "--proteinList", help=" Protein Accession eg. XP_ or WP_047256880.1 in a text file separated by newline. ")
parser.add_argument("-l", "--localGenomeList", help=" Genome File name and Protein Accession ")
parser.add_argument("-ld", "--localGenomeDirectory", help=" Path for Local Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-r", "--redundant", help=" To search all assembly type -r A or -r a but for selected number of assembly eg.,5 for each query use -r 5. ")
parser.add_argument("-e", "--ethreshold", help=" E value threshold. Default = 1e-10 ")
parser.add_argument("-n", "--number", help=" Number of Jackhmmer iterations. Default = 3")
parser.add_argument("-g", "--gene", help=" Number of genes for looking up or downstream. Default = 4 ")
parser.add_argument("-t", "--tree", action="store_true", help=" If you want to see flanking genes along with phylogenetic tree, requires ETE3 installation. By default it will not produce. ")
parser.add_argument("-ts", "--tshape", help=" Size of triangle shapes that represent flanking genes, this option only works when -t is used. Default = 12 ")
parser.add_argument("-tf", "--tfontsize", help=" Size of font inside triangles that represent flanking genes, this option only works when -t is used. Default = 4 ")
parser.add_argument("-to", "--tree_order", action="store_true", help=" Generate Output with Tree, and then use the tree order to generate other view. ")
parser.add_argument("-u", "--user_email", required=True, action="append", metavar="RECIPIENT",default=[], dest="recipients", help=" User Email Address (at least one required) ")
parser.add_argument("-api", "--api_key", help=" NCBI API Key, To get this key kindly check https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ ")
parser.add_argument("-o", "--out_prefix", required= True, help=" Any Keyword to define your output eg. MyQuery ")
parser.add_argument("-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
parser.add_argument("-k", "--keep", action="store_true", help=" If you want to keep the intermediate files eg. gff3 use [-k]. By default it will remove. ")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.2.1')
parser.add_argument("-vb", "--verbose", action="store_true", help=" Use this option to see the work progress for each query as stdout. ")
args = parser.parse_args()
parser.parse_args()




timeout = 10
socket.setdefaulttimeout(timeout)

def checkBioPython(): #Checking Biopython Version
	import Bio
	return (Bio.__version__)

from tkinter.font import Font #Font for postscript to-scale output
myFont12 = Font(family="Helvetica", size=12)
myFont7 = Font(family="Helvetica", size=7)

Entrez.email = args.recipients[0] #User email
Entrez.max_tries = 5
Entrez.sleep_between_tries = 60

if not args.localGenomeList:
	if args.api_key:
		Entrez.api_key = args.api_key #Valid API-key allows 10 queries per seconds, which makes the tool run faster
else:
	if args.api_key:
		print('Since FlaGs will use Local Data api_key is not necessary, Thanks!')
		sys.exit()

#Color generator
def random_color(h=None):
	"""Generates a random color in RGB format."""
	if not h:
		c = random.random()
	d = 0.5
	e = 0.5
	return _hls2hex(c, d, e)

def _hls2hex(c, d, e):
	return '#%02x%02x%02x' %tuple(map(lambda f: int(f*255),colorsys.hls_to_rgb(c, d, e)))

def outliner (item):
	if item =='#ffffff':
		return '#bebebe'
	elif item =='#f2f2f2':
		return '#008000'
	elif item =='#f2f2f3':
		return '#000080'
	else:
		return item

if args.redundant: #Search flanking genes in limited or all available GCFs for each query
	if args.redundant.isdigit():
		if int(args.redundant)==0:
			print('Please use -r option correctly, kindly check help or manual and try again. Thanks!')
			sys.exit()
		else:
			pass
	if not args.redundant.isdigit():
		if args.redundant.lower()=='a':
			pass
		else:
			print('Please use -r option correctly, kindly check help or manual and try again. Thanks!')
			sys.exit()

if args.cpu:
	if int(args.cpu)>0:
		core=int(args.cpu)
	else:
		print('Please use number eg, 1,2...')
		sys.exit()

if args.assemblyList:
	if args.redundant:
		print('"-r" option is only works with "-p", please try again with proper command')
		sys.exit()

if args.localGenomeList:
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

if args.tree_order:
	if not args.tree:
		print("Kindly make sure that you are using -t to make this -to argument working")
		sys.exit()

if args.ethreshold:
	evthresh=args.ethreshold
else:
	evthresh="1e-10"

if args.number:
	iters=args.number
else:
	iters="3"

if args.gene:
	if int(args.gene)>0:
		s= str(int(args.gene)+1)
	else:
		print('Please insert positive values, starting from 1')
		sys.exit()
else:
	s= "5"

if args.localGenomeList:
	if args.localGenomeDirectory:
		if os.path.isdir(args.localGenomeDirectory):
			if args.localGenomeDirectory[-1]=='/':
				localDirIn=args.localGenomeDirectory
				print('Local Data path : ', localDirIn, '\n')
			else:
				localDirIn=args.localGenomeDirectory+'/'
				print('Local Data path : ', localDirIn, '\n')
		else:
			print('No directory Found as : '+ args.localGenomeDirectory)
			sys.exit()
	else:
		localDirIn='./'
else:
	if args.localGenomeDirectory:
		print('Please use -l flag to make -ld flag working')
		sys.exit()

if not args.localGenomeList:
	localDir='./'
else:
	localDir=localDirIn

def checkChar(item): #removing characters
	import re
	items=item.replace('\t','').replace(' ','')
	return re.sub("[a-zA-Z0-9_.]","",items)

queryList=[] #Formatting input as a list

#(queryList)
#eg 1. [['WP_019504790.1', 'GCF_000332195.1'], ['WP_028108719.1', 'GCF_000422645.1']]
#eg 2. [['WP_019504790.1'], ['WP_028108719.1']]

if args.localGenomeList:
	with open (args.localGenomeList, 'r') as gList:
		for line in gList:
			if checkChar(line.rstrip().replace(' ',''))=='':
				Line=line.rstrip().replace(' ','').split('\t')
				if len(Line)<2:
					print('Check Input file, Incorrect Format.')
					sys.exit()
				else:
					newFormat=Line[1]+'\t'+Line[0]
					queryList.append(newFormat.split('\t'))
			else:
				print('The submitted query might include characters not found in NCBI protein accessions eg. > , # , ! etc. Please provide correct format, Thanks!')
				sys.exit()
else:
	if args.proteinList and not args.assemblyList:
		with open (args.proteinList, 'r') as pList:
			for line in pList:
				if checkChar(line.rstrip().replace(' ',''))=='':
					Line=line.rstrip().replace(' ','').split('\t')
					if len(Line)>1:
						print('Check Input file, Incorrect Format.')
						sys.exit()
					else:
						queryList.append(Line)
				else:
					print('The submitted query might include characters not found in NCBI protein accessions eg. > , # , ! etc. Please provide correct format, Thanks!')
					sys.exit()
	elif args.assemblyList and not args.proteinList :
		with open (args.assemblyList, 'r') as apList:
			for line in apList:
				if checkChar(line.rstrip().replace(' ',''))=='':
					Line=line.rstrip().replace(' ','').split('\t')
					if len(Line)<2:
						print('Check Input file, Incorrect Format.')
						sys.exit()
					else:
						newFormat=Line[1]+'\t'+Line[0]
						queryList.append(newFormat.split('\t'))
				else:
					print('The submitted query might include characters not found in NCBI protein accessions eg. > , # , ! etc. Please provide correct format, Thanks!')
					sys.exit()
	else:
		print('Incorrect Input!')
		sys.exit()

def accession_from_xp(accession_nr):

	"""
	:param accession_nr: NCBI protein accession
	:return: Bioproject number of all species for that protein which is used to grab Assembly number
	"""
	#Entrez.email = "_@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
	# blocked, warning is sent to this email.
	try:
		#time.sleep(1)
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
	if bio:
		return set(bio)
	else:
		return {'NAI'}

def accession_from_wp(accession_nr):

	"""
	:param accession_nr: NCBI protein accession
	:return: Set of assembly number of all species for particular protein

	"""
	#Entrez.email = "_@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
	# blocked, warning is sent to this email.
	try:
		#time.sleep(1)
		handle = Entrez.efetch(db="protein", id=accession_nr, rettype="ipg", retmode="xml")
		if float(checkBioPython())<=1.72:
			record = list(Entrez.parse(handle))
			handle.close()
			item=str(record).split(',')
			assembly=[]
			for elements in item:
				if 'assembly' in elements:
					assemblyId=elements.split(':')[1].replace(')','').replace('}','').replace(']','').replace(' ','').replace("'",'')
					assembly.append(assemblyId)
			if assembly and len(assembly)>0:
				return (set(assembly))
			else:
				return {'NAI'}
		else:
			record = Entrez.read(handle)
			handle.close()
			assembly=[]
			item=str(record['IPGReport']).split(',')
			for elements in item:
				if 'assembly' in elements:
					assemblyId=elements.split(':')[1].replace(')','').replace('}','').replace(']','').replace(' ','').replace("'",'')
					assembly.append(assemblyId)
			if assembly and len(assembly)>0:
				return (set(assembly))
			else:
				return {'NAI'}
	except Exception as e:
		print("\t\tQuery {}, not found in database. \n" "\t\tContinuing with the next protein in the list ... \n".format(accession_nr))
		return False

def seq_from_wp(accession_nr):
	"""
	:param accession_nr: NCBI protein accession
	:return: Protein Sequence
	"""
	if accession_nr[-1]!='*':
		#Entrez.email = "_@gmail.com"  # If you do >3 entrez searches on NCBI per second, your ip will be
		# blocked, warning is sent to this email.
		try:
			#time.sleep(1)
			handle = Entrez.efetch(db="protein", id=accession_nr, rettype="gbwithparts", retmode="text")
		except Exception as e:
			print(str(e), ", error in entrez-fetch protein accession, {}, not found in database. \n" "Continuing with the next protein in the list. \nError in function: {}".format(accession_nr, seq_from_wp.__name__))
			return False

		record = SeqIO.read(handle, "genbank")
		handle.close()
		return record.description.split('[')[0]+'\t'+record.seq
	else:
		return accession_nr[:-1]+'\t'+'--'

def remBadChar(item): #removing characters from species name
	import re
	return re.sub("[^a-zA-Z0-9]"," ",item).replace(" ","_")

def identicalProtID(accnr): #searching for identical proteins
	try:
		#time.sleep(1)
		epost_1 = Entrez.read(Entrez.epost(db="protein", id=accnr))
		webenv = epost_1["WebEnv"]
		query_key = epost_1["QueryKey"]
		iden_prots = Entrez.efetch(db="protein", rettype='ipg', retmode='text', webenv=epost_1["WebEnv"], query_key=epost_1["QueryKey"])
		iAccSet=set()
		sAccSet=set()
		for item in iden_prots:
			if item[0:2]!='Id':
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				iAssembly=itemLine[-1]
				if iAccession!=accnr and iAccession[2]=='_':
					iAccSet.add(iAccession)
				if iAccession==accnr and iAccession[2]=='_':
					sAccSet.add(iAccession)
		if len(iAccSet)>0 and len(sAccSet)==0:
			for ids in random.sample(iAccSet,1):
				return ids
		else:
			return accnr
	except:
		return accnr

def identicalProtID_WP(accnr): #searching for identical proteins
	try:
		#time.sleep(1)
		epost_1 = Entrez.read(Entrez.epost(db="protein", id=accnr))
		webenv = epost_1["WebEnv"]
		query_key = epost_1["QueryKey"]
		iden_prots = Entrez.efetch(db="protein", rettype='ipg', retmode='text', webenv=epost_1["WebEnv"], query_key=epost_1["QueryKey"])
		iAccSet=set()
		for item in iden_prots:
			if item[0:2]!='Id':
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				iAssembly=itemLine[-1]
				if iAccession!=accnr and iAccession[:3]=='WP_':
					iAccSet.add(iAccession)
		if len(iAccSet)>0:
			for ids in random.sample(iAccSet,1):
				return ids
		else:
			return accnr
	except:
		return accnr

def identicalProtID_WP_Sp(accnr): #searching for identical proteins with same assembly  for 'NP_417570.1' > WP_000785722.1|GCF_000005845.2
	try:
		#time.sleep(1)
		epost_1 = Entrez.read(Entrez.epost(db="protein", id=accnr))
		webenv = epost_1["WebEnv"]
		query_key = epost_1["QueryKey"]
		iden_prots = Entrez.efetch(db="protein", rettype='ipg', retmode='text', webenv=epost_1["WebEnv"], query_key=epost_1["QueryKey"])
		iAccSetSpecial=set()
		iAssemblyList=[]
		for item in iden_prots:
			if item[0:2]!='Id':
				itemLine=item.rstrip().split('\t')
				iAccession=itemLine[6]
				iAssembly=itemLine[-1]
				if iAccession==accnr and iAccession[2]=='_':
					iAssemblyList.append(iAssembly)
				if iAccession!=accnr and iAccession[:3]=='WP_':
					if iAssemblyList:
						if iAssembly==iAssemblyList[0]:
							assWp=iAccession+'|'+iAssembly
							iAccSetSpecial.add(assWp)
		if iAccSetSpecial:
			for ids in random.sample(iAccSetSpecial,1):
				return ids
		else:
			return '#'
	except:
		return '#'

def sortGCFvsGCA(gcagcfSet):
	if gcagcfSet!='NAI':
		Aset=set()
		Fset=set()
		for items in gcagcfSet:
			if items[2]=='A':
				Aset.add(items)
			if items[2]=='F':
				Fset.add(items)
		if len(Fset)>0:
			return Fset
		elif len(Fset)==0 and len(Aset)>0:
			return Aset
		else:
			return gcagcfSet
	else:
		gcagcfSet

def des_check(item):
	if item:
		return item
	else:
		return 'notFound'

def normalize_strand(item1, item2):  #Strand direction change
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

def lcheck(item):
	if 1 in item:
		return 1
	else:
		return 0

def postscriptSize(item):
	if int(item)<1000:
		return(0)
	else:
		return(int(item)/1000)


def getSpeciesFromGCF(faa):
	import urllib.request
	import requests
	trial=0
	retry=True
	spName=''
	while (retry):
		try:
			#time.sleep(1)
			urlLink="https://www.ncbi.nlm.nih.gov/assembly/"+faa+"/?report=xml&format=text"
			spsname=''
			stsname=''
			fp = urllib.request.urlopen(urlLink)
			for line in fp:
				Line = line.decode("utf8").rstrip()
				if 'SpeciesName' in Line:
					spsname=Line.split(';')[2].split('&')[0]
				if 'Sub_value' in Line:
					stsname=Line.split(';')[2].split('&')[0]
			spName=str(spsname+' '+stsname)
			if spName!=' ':
				retry=False
			else:
				trial+=1
				if args.verbose:
					print('Retrying... For species information of '+ faa + ' # '+ str(trial))
				time.sleep(1)
				if trial % 5 == 0:
					retry=False
				else:
					retry=True
		except:
			trial+=1
			if args.verbose:
				print('Retrying... For species information of '+ faa + ' # '+ str(trial))
			time.sleep(1)
			if trial % 5 == 0:
				retry=False
			else:
				retry=True
	if spName==' ':
		return 'Nothing'
	else:
		if args.redundant:
			return remBadChar(spName)+'_'+remBadChar(faa)
		else:
			return remBadChar(spName)


def spLocal(faa,acc): #getting species name using assembly number or accession
	if faa in speciesNameFromOnlineDict:
		return speciesNameFromOnlineDict[faa]
	else:
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc:
				if args.redundant:
					return remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(faa)
				else:
					return remBadChar(record.description.split('[')[-1][:-1])

def desLocal(faa,acc):
	faaFile=faa+'.faa.gz'
	fastaSeq = gzip.open(localDir+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==acc:
			return record.description.split('[')[0]

def seqLocal(faa,acc):
	faaFile=faa+'.faa.gz'
	fastaSeq = gzip.open(localDir+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==acc:
			return str(record.seq)

def localNone(item):
	if item==None:
		return '--'
	else:
		return item

def seqFasLocal(faa,acc): #making fasta file from accession
	if spLocal(faa,acc):
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc.split('#')[0]:
				record.id=acc+'|'+spLocal(faa,acc)#.replace(':','_').replace('[','_').replace(']','_')
				record.description=''
				return record.format("fasta")
	else:
		faaFile=faa+'.faa.gz'
		fastaSeq = gzip.open(localDir+faaFile, "rt")
		for record in SeqIO.parse(fastaSeq, "fasta"):
			if record.id==acc.split('#')[0]:
				if args.redundant:
					record.id=acc+'|'+remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(faa)#.replace(':','_').replace('[','_').replace(']','_')
				else:
					record.id=acc+'|'+remBadChar(record.description.split('[')[-1][:-1])
				record.description=''
				return record.format("fasta")

def redundantCreate(setDict,nums):
	if nums=='A' or nums=='a':
		newList=random.sample(setDict,len(setDict))
	else:
		if len(setDict)>int(nums):
			newList=random.sample(setDict,int(nums))
		else:
			newList=random.sample(setDict,len(setDict))
	return newList

#Download assembly summary report from NCBI Refseq and genBank
if not args.localGenomeList:
	refDb='./refSeq.db'
	genDb='./genBank.db'
	if os.path.isfile(refDb):
		refDbSize=os.path.getsize(refDb)
	else:
		refDbSize='0'
	if os.path.isfile(genDb):
		genDbSize=os.path.getsize(genDb)
	else:
		genDbSize='0'

	ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
	ftp.cwd("/genomes/refseq") # move to refseq directory

	filenames = ftp.nlst() # get file/directory names within the directory
	if 'assembly_summary_refseq.txt' in filenames:
		ftp.sendcmd("TYPE i")
		if int(ftp.size('assembly_summary_refseq.txt'))!=int(refDbSize):#check if the previously downloaded db exists and if that's updated to recent one
			ftp.retrbinary('RETR ' + 'assembly_summary_refseq.txt', open('refSeq.db', 'wb').write) # get the assembly summary from refseq
		else:
			pass

	ftp_gen = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
	ftp_gen.cwd("/genomes/genbank") # move to refseq directory

	filenames = ftp_gen.nlst() # get file/directory names within the directory
	if 'assembly_summary_genbank.txt' in filenames:
		ftp_gen.sendcmd("TYPE i")
		if int(ftp_gen.size('assembly_summary_genbank.txt'))!=int(genDbSize):#check if the previously downloaded db exists and if that's updated to recent one
			ftp_gen.retrbinary('RETR ' + 'assembly_summary_genbank.txt', open('genBank.db', 'wb').write) # get the assembly summary from refseq
		else:
			pass

	assemblyName={}
	bioDict={} #bioproject as keys and assemble number (eg.GCF_000001765.1) as value
	accnr_list_dict={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
	with open('refSeq.db', 'r') as fileIn:
		for line in fileIn:
			if line[0]!='#':
				Line=line.rstrip().split('\t')
				accnr_list_dict[Line[0]]= Line[7]+'\t'+Line[19]
				bioDict[Line[1]]=Line[0]
				assemblyName[Line[0]]=Line[0]

	ftp_gen = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
	ftp_gen.cwd("/genomes/genbank") # move to refseq directory


	assemblyName_GCA={}
	bioDict_gen={}
	accnr_list_dict_gen={} #create a dictionary accessionNumber is a key and Organism name and ftp Gff3 download Link as value
	with open('genBank.db', 'r') as fileIn:
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

	print ('\n'+ '>> Database Downloaded. Cross-checking of the accession list in progress ...'+ '\n')


q=0
ne=0
queryDict={} #protein Id as query and a set of assembly number as value [either All or Species of interest]
#{'WP_019504790.1#1': {'GCF_000332195.1'}, 'WP_028108719.1#2': {'GCF_000422645.1'}, 'WP_087820443.1#3': {'GCF_900185565.1'}}
#{'WP_019504790.1#1': 'GCF_000332195.1', 'WP_028108719.1#2': 'GCF_000422645.1', 'WP_087820443.1#3': 'GCF_900185565.1'} local
if not args.localGenomeList:
	with open (args.out_prefix+'_NameError.txt', 'w') as fbad:
		for query in queryList:
			q+=1
			accession_from_wp_out=''
			accession_from_xp_out=''
			identicalProtID_Out=''
			accession_from_wp_ID_out=''
			accession_from_xp_ID_out=''
			accession_from_wp_IDSame_out=''
			exceptionalWP_out=''
			accession_from_wp_exceptional=''
			special_out = ''
			if args.verbose:
				print('\t Checking Query '+ query[0] +' ....'+ '('+str(q)+'/'+str(len(queryList))+')')
			if len(query)<2:
				if query[0][:2]=='WP' and query[0][-2]=='.': #WP Accession full WP_000785722.1
					accession_from_wp_out=accession_from_wp(query[0])
					if accession_from_wp_out:
						queryDict[query[0]+'#'+str(q)]=sortGCFvsGCA(accession_from_wp_out)
					else:
						ne+=1
						print(query[0], file= fbad)
				elif query[0][:2]=='XP' and query[0][-2]=='.': #XP Accession full XP_003256407.1
					accession_from_xp_out=accession_from_xp(query[0])
					if accession_from_xp_out:
						assemList=[]
						for bioprojs in accession_from_xp_out:
							if bioprojs in bioDict:
								assemList.append(bioDict[bioprojs])
						if assemList:
							queryDict[query[0]+'#'+str(q)]=sortGCFvsGCA(set(assemList))
					else:
						ne+=1
						print(query[0], file= fbad)
				else: #other accessions can be XP_003256407 , WP_000785722 too
					identicalProtID_Out=identicalProtID(query[0])
					if identicalProtID_Out!=query[0]: #can be anything XP_ or WP_ or YP_ or NP_
						if identicalProtID_Out[:-3]!='XP_': #Not XPs
							accession_from_wp_ID_out=accession_from_wp(identicalProtID_Out)
							if accession_from_wp_ID_out:
								asset=set()
								for elements in accession_from_wp_ID_out:
									asset.add(elements)
								if len(asset)>0:
									queryDict[identicalProtID_Out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
							else:
								ne+=1
								print(query[0], file= fbad)
						if identicalProtID_Out[:-3]=='XP_': #if XPs
							accession_from_xp_ID_out=accession_from_xp(identicalProtID_Out)
							if accession_from_xp_ID_out:
								assemList=[]
								for bioprojs in accession_from_xp_ID_out:
									if bioprojs in bioDict:
										assemList.append(bioDict[bioprojs])
								if assemList:
									queryDict[identicalProtID_Out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(set(assemList))
							else:
								ne+=1
								print(query[0], file= fbad)
					elif identicalProtID_Out==query[0]: #excluding XP Wp pre  #YP NP
						exceptionalWP_out = identicalProtID_WP(identicalProtID_Out)
						special_out = identicalProtID_WP_Sp(identicalProtID_Out) #list Query GCF
						if not args.redundant:
							if special_out!='#':
								asset=set()
								asset.add(special_out.split('|')[1])
								queryDict[query[0]+'#'+str(q)]=asset
							else:
								if exceptionalWP_out[:3]=='WP_':
									accession_from_wp_exceptional=accession_from_wp(exceptionalWP_out)
									if accession_from_wp_exceptional:
										asset=set()
										for elements in accession_from_wp_exceptional:
											asset.add(elements)
										if len(asset)>0:
											queryDict[exceptionalWP_out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
									else:
										ne+=1
										print(query[0], file= fbad)
						if args.redundant:
							if exceptionalWP_out[:3]=='WP_':
								accession_from_wp_exceptional=accession_from_wp(exceptionalWP_out)
								if accession_from_wp_exceptional:
									asset=set()
									for elements in accession_from_wp_exceptional:
										asset.add(elements)
									if len(asset)>0:
										queryDict[exceptionalWP_out+'#'+str(q)+'.'+query[0]]=sortGCFvsGCA(asset)
								else:
									ne+=1
									print(query[0], file= fbad)
			else:
				if query[0][:3]=='XP_' and query[0][-2]=='.': #XP Accession
					asset=set()
					accession_from_xp_out=accession_from_xp(query[0])
					if accession_from_xp_out:
						for bioprojs in accession_from_xp_out:
							if bioprojs in bioDict:
								if bioDict[bioprojs]==query[1]:
									asset.add(query[1])
						if len(asset)>0:
							queryDict[query[0]+'#'+str(q)]=asset
					else:
						ne+=1
						print(query[0], file= fbad)
				elif query[0][:3]!='XP_' and query[0][-2]=='.': #not XP Accession
					asset=set()
					accession_from_wp_out=accession_from_wp(query[0])
					if accession_from_wp_out:
						for elements in accession_from_wp_out:
							if query[1]==elements:
								asset.add(query[1])
					if len(asset)>0:
						queryDict[query[0]+'#'+str(q)]=asset
				else:
					ne+=1
					print(query[0], file= fbad)
else:
	for query in queryList:
		q+=1
		queryDict[query[0]+'#'+str(q)]=query[1]




nai=0
NqueryDict={} #{'WP_019504790.1#1': ['GCF_000332195.1'], 'WP_028108719.1#2': ['GCF_000422645.1'], 'WP_087820443.1#3': ['GCF_900185565.1']}
if args.localGenomeList:
	with open (args.out_prefix+'_Insufficient_Info_In_DB.txt', 'w') as fNai:
		for query in queryDict:
			assemblyIdlist=[]
			if queryDict[query]!={'NAI'}:
				assemblyId=queryDict[query]
				faa_gz=localDir+queryDict[query]+'.faa.gz'
				if os.path.isfile(faa_gz):
					with gzip.open(localDir+assemblyId+'.faa.gz', 'rb') as faaIn:
						for line in faaIn:
							if line.decode('utf-8')[0]=='>':
								Line=line.decode('utf-8').rstrip()
								if '>'+query.split('#')[0]==Line.split(' ')[0]:
									gff_gz=localDir+assemblyId+'.gff.gz'
									if os.path.isfile(gff_gz):
										with gzip.open(localDir+assemblyId+'.gff.gz', 'rb') as lgffIn: #Download and read gff.gz
											cds_c=0
											name_c=0
											prot_c=0
											cset=set()
											nset=set()
											pset=set()
											for line in lgffIn:
												if line.decode('utf-8')[0]!='#':
													Line=line.decode('utf-8').rstrip().split('\t')
													if Line[2]=='CDS':
														cds_c=1
														cset.add(cds_c)
														if Line[8].split(';')[3][:5]=='Name=': #eliminates pseudo gene as they don't have 'Name='
															name_c=1
															nset.add(name_c)
															if Line[8].split(';')[3].split('=')[1]==query.split('#')[0]:
																assemblyIdlist.append(assemblyId)
																NqueryDict[query]=list(set(assemblyIdlist))
																prot_c=1
																pset.add(prot_c)
															else:
																prot_c=0
																pset.add(prot_c)
														else:
															name_c=0
															nset.add(name_c)
													else:
														cds_c=0
														cset.add(cds_c)
											if lcheck(pset)>0:
												pass
											elif lcheck(pset)==0:
												if lcheck(cset)>0 and lcheck(nset)>0:
													print(query.split('#')[0],' did not match with supplement local GFF File')
													print(query, file=fNai)
													nai+=1
											else:
												print('Use recommended [NCBI refseq] format of GFF file')
												break
									else:
										print("Error: %s file not found" % gff_gz)
				else:
					print("Error: %s file not found" % faa_gz)
			else:
				print(query, file=fNai)
				nai+=1
else:
	with open (args.out_prefix+'_Insufficient_Info_In_DB.txt', 'w') as fNai:
		for query in queryDict:
			#print('\t', query, queryDict[query], 'All')
			if queryDict[query]:
				if len(queryDict[query])!=0:
					if queryDict[query]!={'NAI'}:
						#print('\t', query, queryDict[query], 'filtered')
						if args.redundant:
							redun=0
							for newRed in (redundantCreate(queryDict[query],args.redundant)):
								redun+=1
								NqueryDict[query+'.'+str(redun)]=list(str(newRed).split())
						else:
							NqueryDict[query]=random.sample(queryDict[query],1)
					else:
						print(query, file=fNai)
						nai+=1
				else:
					print(query, file=fNai)
					nai+=1
			else:
				print(query, file=fNai)
				nai+=1

if not args.localGenomeList:
	print('\n> Downloading Genome Assembly Files from NCBI FTP Server \n')




newQ=0
speciesNameFromOnlineDict={}
for query in NqueryDict:
	newQ+=1
	total=0
	for item in NqueryDict[query]:
		speciesNameFromOnlineDict[item]=getSpeciesFromGCF(item)
		if not args.localGenomeList:
			AssemDown=0
			AssemFailed=0
			if item in accnr_list_dict:
				retry = True
				while (retry):
					try:
						ftpLine=accnr_list_dict[item].split('\t')[1]
						ftpsplitDir = ftpLine.split('/')[3:]
						ftp_path = '/'.join(map(str,ftpsplitDir))
						#time.sleep(1)
						ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
						ftp.cwd('/'+ftp_path)
						files = ftp.nlst()
						FileToDownload=[]
						for ftpelements in files:
							if '_genomic.gff.gz' in ftpelements:
								FileToDownload.append(ftpelements)
							if '_protein.faa.gz' in ftpelements:
								FileToDownload.append(ftpelements)
						if len(FileToDownload)==2:
							for elements in FileToDownload:
								#print(elements) #GCF_000964005.1_WiARP1.0_genomic.gff.gz GCF_000964005.1_WiARP1.0_protein.faa.gz
								if '_genomic.gff.gz' in elements: #Check if GFF.gz is there
									if args.verbose:
										ftp.set_debuglevel(1)
									ftp.set_pasv(True)
									ftp.voidcmd('TYPE I')
									ftp.sendcmd("TYPE i")
									try:
										gffFileSize=ftp.size(elements)
										gffFileName=item+'.gff.gz'
										gffFileDownloaded=open(item+'.gff.gz', 'wb')
										ftp.retrbinary('RETR ' + elements, gffFileDownloaded.write)
										ftp.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
										ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 15)
										ftp.sock.setsockopt(socket.IPPROTO_TCP,socket.TCP_KEEPCNT, 8)
										ftp.voidcmd("NOOP")
										gffFileDownloaded.close()
										if os.path.isfile(localDir+gffFileName):
											if os.path.getsize(gffFileName)==gffFileSize:
												AssemDown+=1
												retry = False
									except:
										AssemFailed+=-1
										retry = True
								if '_protein.faa.gz' in elements: #Check if faa.gz is there
									if args.verbose:
										ftp.set_debuglevel(1)
									ftp.set_pasv(True)
									ftp.voidcmd('TYPE I')
									ftp.sendcmd("TYPE i")
									ftp.sendcmd("TYPE i")
									try:
										faaFileSize=ftp.size(elements)
										faaFileName=item+'.faa.gz'
										faaFileDownloaded=open(item+'.faa.gz', 'wb')
										ftp.retrbinary('RETR ' + elements, faaFileDownloaded.write)
										ftp.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
										ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 15)
										ftp.sock.setsockopt(socket.IPPROTO_TCP,socket.TCP_KEEPCNT, 8)
										ftp.voidcmd("NOOP")
										faaFileDownloaded.close()
										if os.path.isfile(localDir+faaFileName):
											if os.path.getsize(faaFileName)==faaFileSize:
												AssemDown+=1
												retry = False
									except:
										AssemFailed+=-1
										retry = True
						else:
							AssemFailed+=-2
							retry = False
					except:
						retry = True
			else:
				AssemFailed+=-2
			total=AssemDown+AssemFailed
			if total==2:
				if args.verbose:
					print('\n\t'+'> NCBI Genome Assembly has been downloaded for '+query.split('#')[0]+' ('+str(newQ)+'/'+str(len(NqueryDict))+')'+'\n')

for query in NqueryDict:
	for item in NqueryDict[query]:
		if speciesNameFromOnlineDict[item]=='Nothing':
			faaFile=item+'.faa.gz'
			fastaSeq = gzip.open(localDir+faaFile, "rt")
			for record in SeqIO.parse(fastaSeq, "fasta"):
				if record.id==query.split('#')[0]:
					if args.redundant:
						speciesNameFromOnlineDict[item]=remBadChar(record.description.split('[')[-1][:-1])+'_'+remBadChar(item)
					else:
						speciesNameFromOnlineDict[item]=remBadChar(record.description.split('[')[-1][:-1])


if args.keep:
	with open(args.out_prefix+'_speciesInfo.txt','w') as asmOut:
		for query in NqueryDict:
			for item in NqueryDict[query]:
				print(item, query.split('#')[0], speciesNameFromOnlineDict[item], sep='\t', file=asmOut)

#print(NqueryDict) #{'WP_019504790.1#1': ['GCF_000332195.1'], 'WP_028108719.1#2': ['GCF_000422645.1'], 'WP_087820443.1#3': ['GCF_900185565.1']}

print('\n'+'>> Input file assessment report: ')
print('\t'+'Discarded protein ids with improper accession : '+str(ne)+'. See "'+args.out_prefix+'_NameError.txt'+'" file for details.')
print('\t'+'Discarded protein ids lacking proper information in RefSeq DB : '+str(nai)+'. See "'+args.out_prefix+'_Insufficient_Info_In_DB.txt'+'" file for details.')
print('\t'+'Remaining queries: '+str(newQ))



def getGeneId(item):
	matchObject = re.search('(GeneID:.*?,)', item)
	if matchObject:
		return matchObject.group(1)[:-1]

def getGeneId_gene(item):
	matchObject2 = re.search('(GeneID:.*;)', item)
	if matchObject2:
		return matchObject2.group(1).split(';')[0]


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
treeFastadict={} #Query as key and sequence in fasta as value
querySeqDict={} #For Tree Command


count=0
for query in NqueryDict:
	count+=1
	if args.verbose:
		print('\n'+'> '+str(count)+' in process out of '+str(newQ)+' ... '+'\n')
		print('Query Name =', query.split('#')[0], '\n')
	for item in NqueryDict[query]:  #{'WP_019504790.1#1': ['GCF_000332195.1'],NqueryDict
		a=0
		LineList=[]
		geneProt={} # 'gene2504': 'WP_092248795.1', 'gene1943': 'tRNA'
		geneChrom={} #'gene1708': 'NZ_MJLP01000030.1'
		#item='GCF_'+items[4:]
		speciesNameFromDB=speciesNameFromOnlineDict[item]
		gff_gz=localDir+item+'.gff.gz'
		if os.path.isfile(gff_gz):
			LineList=[]
			geneProt={} # 'gene2504': 'WP_092248795.1', 'gene1943': 'tRNA'
			geneChrom={} #'gene1708': 'NZ_MJLP01000030.1'
			with gzip.open(localDir+item+'.gff.gz', 'rb') as gffIn: #Download and read gff.gz
				for line in gffIn:
					if line.decode('utf-8')[0]!='#':
						Line=line.decode('utf-8').rstrip().split('\t')
						if Line[2]=='CDS':
							if Line[8].split(';')[3][:5]=='Name=': #eliminates pseudo gene as they don't have 'Name='
								if query.split('_')[0]=='XP':
									if 'GeneID:' in Line[8]:
										geneProt[getGeneId(Line[8])]=Line[8].split(';')[3].split('=')[1]
										geneChrom[getGeneId(Line[8])]=Line[0]
										#print(getGeneId(Line[8]), Line[8].split(';')[3].split('=')[1], Line[0], 'gpc#')
								else:
									#print(Line[8].split(';')[1].split('=')[1], Line[8].split(';')[3].split('=')[1], Line[0], 'gpc#')
									geneProt[Line[8].split(';')[1].split('=')[1]]=Line[8].split(';')[3].split('=')[1]
									geneChrom[Line[8].split(';')[1].split('=')[1]]=Line[0]
									##print(geneProt)>'GeneID:187667': 'NP_493855.2' NP_417570.1 gene-b3099
						if Line[2][-4:]=='gene':
							a+=1
							if query.split('_')[0]=='XP':
								if 'GeneID:' in Line[8]:
									newGene=str(a)+'\t'+getGeneId_gene(Line[8])+'\t'+ Line[3]+'\t'+Line[4]+'\t'+ Line[6]+ '\t'+ Line[0]
									#print(newGene) #1	GeneID:353377	3747	3909	-	NC_003279.8
									LineList.append(newGene.split('\t'))
									for genDes in Line[8].split(';'):
										if 'gene_biotype=' in genDes:
											if getGeneId_gene(Line[8]) not in geneProt:
												geneProt[getGeneId_gene(Line[8])]=genDes.split('=')[1]+'_'+query.split('#')[1]+'.'+str(random.randint(0,int(s)*2-1))+'*'
							else:
								newGene=str(a)+'\t'+Line[8].split(';')[0][3:]+'\t'+ Line[3]+'\t'+Line[4]+'\t'+ Line[6]+ '\t'+ Line[0]
								LineList.append(newGene.split('\t')) #1       gene3006        10266   10342   -       NZ_FPCC01000034.1
								for genDes in Line[8].split(';'):
									if 'gene_biotype=' in genDes:
										if Line[8].split(';')[0][3:] not in geneProt:
											geneProt[Line[8].split(';')[0][3:]]=genDes.split('=')[1]+'_'+query.split('#')[1]+'.'+str(random.randint(0,int(s)*2-1))+'*'
				geneList=[] ##List of gene names coding same protein (accession), we are taking one from them
				for genes in geneProt:
					if geneProt[genes]==query.split('#')[0]:
						geneList.append(genes)
				if len(geneList)>0:
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
										if args.tree:
											treeFastadict[query]=str(seqFasLocal(item,query))
											querySeqDict[query+'|'+remBadChar(spLocal(item,query.split('#')[0]))]=str(seqLocal(item,query.split('#')[0]))
										if speciesNameFromDB!='Nothing' or speciesNameFromDB!='':
											#print(speciesNameFromDB, 1)
											speciesDict[query]=speciesNameFromDB
										else:
											#print(spLocal(item, query.split('#')[0]), 2)
											speciesDict[query]=spLocal(item, query.split('#')[0])
										queryStrand[query]= LineList[LineList.index(line)][4]
										positionDict[query]= ("\t".join(map(str,LineList[LineList.index(line)][2:-2])))
										LengthDict[query]= int(LineList[LineList.index(line)][3])-int(LineList[LineList.index(line)][2])+1
										udsDict={}
										dsDict={}
										udsDict[0]= query+'+' # O strand
										for x in range(1,int(s)):
											if LineList.index(line)-x>=0 and LineList.index(line)-x<len(LineList):
												if rangeList.count(int(LineList[LineList.index(line)-x][0]))>0:
													acc_CGF_Dict[query]= LineList[LineList.index(line)-x][-1] +'\t'+ item
													seqDict[str(geneProt[LineList[LineList.index(line)-x][1]])]=localNone(seqLocal(item, geneProt[LineList[LineList.index(line)-x][1]]))
													desDict[geneProt[LineList[LineList.index(line)-x][1]]]=desLocal(item, geneProt[LineList[LineList.index(line)-x][1]])
													positionDict[geneProt[LineList[LineList.index(line)-x][1]]+'#'+query.split('#')[1]]= ("\t".join(map(str,LineList[LineList.index(line)-x][2:-2])))
													LengthDict[geneProt[LineList[LineList.index(line)-x][1]]+'#'+query.split('#')[1]]= int(LineList[LineList.index(line)-x][3])-int(LineList[LineList.index(line)-x][2])+1
													udsDict[int(ups(LineList[LineList.index(line)][4])[0]+str(x))]= geneProt[LineList[LineList.index(line)-x][1]]+'#'+query.split('#')[1]+\
														normalize_strand(LineList[LineList.index(line)][4],LineList[LineList.index(line)-x][4])
										for y in range(1,int(s)):
											if LineList.index(line)+y<len(LineList):
												if rangeList.count(int(LineList[LineList.index(line)+y][0]))>0:
													acc_CGF_Dict[query]= LineList[LineList.index(line)+y][-1] +'\t'+ item
													seqDict[str(geneProt[LineList[LineList.index(line)+y][1]])]=localNone(seqLocal(item,geneProt[LineList[LineList.index(line)+y][1]]))
													desDict[geneProt[LineList[LineList.index(line)+y][1]]]=desLocal(item,geneProt[LineList[LineList.index(line)+y][1]])
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
				else:
					FlankFoundDict[query]='No'
					FoundDict[query]='No: ProteinID was not found in Genome Assembly'
					if args.verbose:
						print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')
		else:
			FlankFoundDict[query]='No'
			FoundDict[query]='No: ProteinID was not found in Genome Assembly'
			if args.verbose:
				print('\t', query.split('#')[0], 'Report: Flanking Genes Not Found', '\n')

if not args.localGenomeList:
	if args.keep:
		pass
	else:
		subprocess.Popen("rm GC*_*.gz", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

allFlankGeneList=[]
for keys in accFlankDict:
	for item in accFlankDict[keys]:
		allFlankGeneList.append(accFlankDict[keys][item].split('#')[0])


flankF=0
with open (args.out_prefix+'_flankgene_Report.log', 'w') as errOut:
	serial=0
	print('#Serial','Query','Assembly_Found', 'FlankingGene_Found', sep='\t', file=errOut)
	for queries in NqueryDict:
		serial+=1
		if queries in FoundDict:
			if queries in FlankFoundDict:
				if FlankFoundDict[queries]=='Yes':
					flankF+=1
				print(str(serial), queries.split('#')[0], FoundDict[queries], FlankFoundDict[queries], sep='\t', file=errOut)
			else:
				print(str(serial), queries.split('#')[0], FoundDict[queries], 'No', sep='\t', file=errOut)
		else:
			print(str(serial), queries.split('#')[0], 'No', 'No', sep='\t', file=errOut)

reportDict={}
for query in queryList:
	queryNumList=[]
	for queryNum in NqueryDict:
		if query[0] in queryNum:
			queryNumList.append(queryNum)
	if len(queryNumList)>0:
		reportDict[query[0]]=queryNumList
	else:
		reportDict[query[0]]=str('no').split()

def similarityID (item1, item2):
	if item1==item2:
		return 'Same'
	else:
		return 'Changed'

def reporter(i1,i2,i3,i4,i5):
	if i2=='No' and i3=='No' and i4=='No' and i5=='No':
		return i1 + ' Failed :  No record for accession ' + i1
	if i2!='No' and i3=='Same' and i4!='No' and i5=='No':
		return i1 + ' is a valid NCBI protein accession but Discarded :  No Flanking Gene was found for ' + i1 + ' in Assembly ID '+ i4
	if i2!='No' and i3=='Same' and i4!='No' and i5!='No':
		return i1 + ' is valid as a NCBI protein accession and reported in Assembly ID '+ i4
	if i2!='No' and i3=='Changed' and i4!='No' and i5=='No':
		return i1 + ' is invalid NCBI protein accession therefore converted to identical RefSeq sequence with accession '+ i2 + ' but Discarded :  No Flanking Gene was found for ' + i2 + 'in Assembly ID '+ i4
	if i2!='No' and i3=='Changed' and i4!='No' and i5!='No':
		return i1 + ' is invalid NCBI protein accession therefore converted to identical RefSeq sequence with accession '+ i2 + ' which is reported in Assembly ID '+ i4

qcount=0
with open (args.out_prefix+'_QueryStatus.txt', 'w') as sumOut:
	print('#Serial', 'Status', sep='\t', file=sumOut)
	for query in queryList:
		qcount+=1
		for item in reportDict[query[0]]:
			if item in FlankFoundDict:
				if FlankFoundDict[item]=='Yes':
					print(str(qcount), reporter(query[0], item.split('#')[0], similarityID(query[0], item.split('#')[0]), ''.join(NqueryDict[item]), 'Yes'), sep='\t', file=sumOut)
				else:
					print(str(qcount), reporter(query[0], item.split('#')[0], similarityID(query[0], item.split('#')[0]), ''.join(NqueryDict[item]), 'No'), sep='\t', file=sumOut)
			else:
				print(str(qcount), reporter(query[0], 'No', 'No', 'No', 'No'), sep='\t', file=sumOut)

print('\n'+'>> Flanking Genes found : '+str(flankF)+' out of remaining '+str(serial)+'. See "'+args.out_prefix+'_flankgene_Report.log'+'" file for details.'+'\n'+'\n')

if int(flankF)==0:
	print('>> No Flanking Genes found, please update your accession list')
	sys.exit()
else:
	pass


if args.tree: #Generate fasta file for making phylogenetic Tree
	with open(args.out_prefix+'_tree.fasta', 'w') as treeOut:
		for queries in NqueryDict:
			if queries in treeFastadict:
				print(treeFastadict[queries], file=treeOut)


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

print('\n>> Now running Jackhmmer and clustering flanking genes\n')

infilename=args.out_prefix+'_flankgene.fasta'+'_cluster_out'

directory = args.out_prefix+'_flankgene.fasta'+'_cluster_out_individuals'
if not os.path.exists(directory):
	os.makedirs(directory)

infile=open(args.out_prefix+'_flankgene.fasta'+'_cluster_out',"r").read()
al=infilename+"_"+iters+"_"+evthresh+"_jackhits.tsv"
outacclists=open(al,"w")

percentileJack=0
i=1
for seqids in sorted(seqDict): #Running Jackhmmer for finding homologs
	if seqDict[seqids]!='--':
		percentileJack+=1
		if args.verbose:
			if percentileJack % 5 == 0:
				print('\t'+'>>> '+str(round(int(percentileJack)*100/b))+'%'+' Completed...'+'('+str(percentileJack)+'/'+str(b)+')')
			if percentileJack % b == 0:
				print('\t'+'>>> Completed ' +'\n')
		i_f=directory+"/"+str(i)+".txt"
		indivfile=open(i_f,"w")
		indivfile.write(">"+seqids+'\n'+seqDict[seqids])
		indivfile.close()
		if args.cpu:
			command="jackhmmer --cpu %s -N %s --incE %s --incdomE %s --tblout %s/tblout%s.txt %s  %s>%s/out%s.txt" %(core, iters, evthresh, evthresh, directory, str(i), i_f, infilename, directory, str(i))
		else:
			command="jackhmmer -N %s --incE %s --incdomE %s --tblout %s/tblout%s.txt %s  %s>%s/out%s.txt" %(iters, evthresh, evthresh, directory, str(i), i_f, infilename, directory, str(i))
		os.system(command)
		tbl=open(directory+"/tblout"+str(i)+".txt").read()
		part=tbl.split("----------\n")[1].split("\n#")[0]
		lines=part.splitlines()
		acclist=[]
		for line in lines:
			lineList=line.split()
			if len(lineList)>17:
				inc=line.split()[17]
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

#d : {1: ['WP_000153877.1'], 2: [], 3: ['WP_000291520.1'], 4: ['WP_000342211.1'], 5: [],... 30: ['WP_055032751.1', 'WP_001246052.1']}

trueAccessionCount={}
for keys in d:
	numbers=[]
	for item in d[keys]:
		numbers.append(allFlankGeneList.count(item))
	trueAccessionCount[(';'.join(map(str,d[keys])))]=sum(numbers)

#trueAccessionCount: 'WP_001229260.1;WP_001229255.1': 4
odtrue=OrderedDict(sorted(trueAccessionCount.items(), key= lambda item:item[1],reverse=True))

familyNumber=0
with open(infilename+"_"+iters+"_"+evthresh+"_clusters.tsv","w") as clusOut:
	for k, v in odtrue.items():
		#print (k.split(';'), odtrue[k], v)
		if len(k.split(';'))>0 and v>0:
			familyNumber+=1
			print(str(familyNumber),str(odtrue[k]),k, sep='\t', file=clusOut)

outfile_des=open(infilename+"_"+iters+"_"+evthresh+"_outdesc.txt","w")
inf=open(infilename+"_"+iters+"_"+evthresh+"_clusters.tsv","r")

acclists=inf.read().splitlines()
for line in acclists:
	acclist=line.split("\t")[2].split(";")
	familyAssignedValue=line.split("\t")[0]
	if int(line.split("\t")[1])>1:
		for acc in acclist:
			outfile_des.write(familyAssignedValue+'('+str(allFlankGeneList.count(acc))+')'+"\t"+acc+"\t"+desDict[acc]+"\n")
		outfile_des.write ("\n\n")

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
noProtP=int(max(familynum))+3
noColor=int(max(familynum))+4
for ids in LengthDict:
	if ids.split('#')[0][-1]=='*':
		if ids.split('#')[0][:2].lower()=='ps':
			familyDict[ids.split('#')[0]]=noProtP
		else:
			familyDict[ids.split('#')[0]]=noProt
	if ids.split('#')[0] not in familyDict:
		if ids in NqueryDict:
			familyDict[ids.split('#')[0]]=center
		else:
			familyDict[ids.split('#')[0]]=noColor

color={}
color[noColor]='#ffffff'
color[center]='#000000'
color[noProt]='#f2f2f2'
color[noProtP]='#f2f2f3'
colorDict={} #Assigned family Number from Jackhammer : colorcode
for families in set(familynum):
	if families == 0:
		colorDict[families]=str('#ffffff')
	else:
		if random_color()!='#ffffff' or random_color()!='#000000' or random_color()!='#f2f2f2' or random_color()!='#f2f2f3' :
			colorDict[families]=random_color()

colorDict.update(color)

maxs=(int(s)-1) # required to calculate border size of postscript output
mins=maxs-(maxs*2) # required to calculate border size of postscript output

if not args.tree_order:
	nPos=[]
	pPos=[]
	with open(args.out_prefix+'_operon.tsv', 'w') as opOut:
		for queries in accFlankDict:
			for items in sorted(accFlankDict[queries]):
				if queryStrand[queries]=='+':
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'|'+remBadChar(speciesDict[queries])
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
					lengths=LengthDict[accFlankDict[queries][items][:-1]] #c2
					species=queries+'|'+remBadChar(speciesDict[queries])
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

	windowMost=round(((max(pPos)+abs(min(nPos))+1)*4)/100)
	widthM=(windowMost*3)+500
	heightM=int(newQ)*20
	from tkinter import *
	master = Tk()

	canvas = Canvas(master, width=widthM,height=heightM,background='white', scrollregion=(0,0,round(widthM*2.5),round(heightM*2.5)))
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
		elif item==noProtP:
			return ' '
		elif item==noColor:
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
			ptnstats=entries[0].split("\t") #c2
			org=ptnstats[0][:ptnstats[0].index('|')]+ptnstats[0][ptnstats[0].index('|'):].replace('_',' ')
			textspace=widthM/2
			line_pos_y=line_pos_y+16-round(postscriptSize(newQ))
			half_dom_height=5-round(postscriptSize(newQ))
			text = canvas.create_text(textspace/2-textspace/8,line_pos_y, text=org, fill="#404040", font=myFont12)
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
				textd1 = canvas.create_text(dom1_start+(dom1_len/2),line_pos_y, text=operonFamily(dom1_name), font=myFont7)
				coln=coln+1

	retval = canvas.postscript(file=args.out_prefix+"_flankgenes.ps", height=heightM, width=widthM, colormode="color")



if args.tree:###Tree Command with ETE###
	tree_file= args.out_prefix+'_tree.fasta'
	if args.cpu:
		tree_command="ete3 build -a %s -o %s --nochecks --clearall -w mafft_default-trimal01-none-fasttree_full --rename-dup-seqnames --cpu %s" %(tree_file, tree_file[:-6], core)
	else:
		tree_command="ete3 build -a %s -o %s --nochecks --clearall -w mafft_default-trimal01-none-fasttree_full --rename-dup-seqnames" %(tree_file, tree_file[:-6])
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
		elif item==noProtP:
			return ' '
		elif item==noColor:
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

	nwTree=''
	motifDict={}
	motifDict_2={}
	if os.path.isfile(args.out_prefix+'_tree/mafft_default-trimal01-none-fasttree_full/'+args.out_prefix+'_tree.fasta.final_tree.nw') == True:
		with open(args.out_prefix+'_tree/mafft_default-trimal01-none-fasttree_full/'+args.out_prefix+'_tree.fasta.final_tree.nw', 'r') as treeIn:
			for line in treeIn:
				nwTree=line
				for items in line.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
					item=items.split('|')[0]
					simple_motifs=[]
					simple_motifs_2=[]
					for keys in sorted(startDict):
						if keys in accFlankDict[item]:
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
	else:
		print('> ETE3 failed to create tree due to lack of valid protein accesions, at least 2 required !')
		sys.exit()

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
			t.render(args.out_prefix+"_flankgenes_2.svg", tree_style=ts)

if args.tree and args.tree_order:  # Queries in postscript file will be presented as tree order
	treeOrderList=[]
	with open(args.out_prefix+'_ladderTree.nw', 'r') as laddertreeIn:
		for line in laddertreeIn:
			for items in line.replace('(','').replace(')', '').replace(';', '').replace(',','\t').split('\t'):
				item=items.split('|')[0]
				treeOrderList.append(item)

	ntPos=[]
	ptPos=[]
	with open(args.out_prefix+'_TreeOrder_operon.tsv', 'w') as opOut:
		for queries in treeOrderList:
			for items in sorted(accFlankDict[queries]):
				if queryStrand[queries]=='+':
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'|'+remBadChar(speciesDict[queries])
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
					ntPos.append(nP)
					ptPos.append(pP)

				else:
					ids=accFlankDict[queries][items][:-1]
					lengths=LengthDict[accFlankDict[queries][items][:-1]]
					species=queries+'|'+remBadChar(speciesDict[queries])
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
					ntPos.append(nP)
					ptPos.append(pP)

			print('\n\n', file=opOut)

	windowMost=round(((max(ptPos)+abs(min(ntPos))+1)*4)/100)
	widthM=(windowMost*3)+500
	heightM=int(newQ)*20
	from tkinter import *
	master = Tk()


	canvas = Canvas(master, width=widthM,height=heightM,background='white', scrollregion=(0,0,round(widthM*2.5),round(heightM*2.5)))
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
		elif item==noProtP:
			return ' '
		elif item==noColor:
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
			org=ptnstats[0][:ptnstats[0].index('|')]+ptnstats[0][ptnstats[0].index('|'):].replace('_',' ')
			textspace=widthM/2
			line_pos_y=line_pos_y+16-round(postscriptSize(newQ))
			half_dom_height=5-round(postscriptSize(newQ))
			text = canvas.create_text(textspace/2-textspace/8,line_pos_y, text=org, fill="#404040", font=myFont12)
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
				textd1 = canvas.create_text(dom1_start+(dom1_len/2),line_pos_y, text=operonFamily(dom1_name), font=myFont7)
				coln=coln+1

	retval = canvas.postscript(file=args.out_prefix+"_treeOrder_flankgenes.ps", height=heightM, width=widthM, colormode="color")

print('\n'+'<<< Done >>>')
print('For FlaGs Citation: https://www.biorxiv.org/content/10.1101/362095v2')
sys.exit()
