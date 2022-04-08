import os
import gzip
import re
import argparse

class Sample:

	def __init__(self, name):
		self.reads = []
		self.name = name
		self.type = ""
				
	def getHisatString(self):
		pass
			
	def addRead(self, read):
		pass
		
	def getName(self):
		return self.name
	
	def __str__(self):
		return self.name + " (" +str(len(self.reads)) + " files - "+ self.type +")"
	



class PairedReadSample(Sample):
	def __init__(self, name):
		super().__init__(name)
		self.left_reads = []
		self.right_reads = []
		self.type = "paired"
		
		
	def getHisatString(self):
		hs = "-1 "
		for i in range(len(self.left_reads)):
			read = self.left_reads[i]
			hs += str(read)
			if i < len(self.left_reads)-1:
				hs += ","	
		
		hs += " -2 "
		for i in range(len(self.right_reads)):
			read = self.right_reads[i]
			hs += str(read)
			if i < len(self.right_reads)-1:
				hs += ","	
		return hs
	
	
		
	def addRead(self, read):
		self.reads.append(read)
		if "R1" in read or "READ1" in read.upper():
			self.left_reads.append(read)
		elif "R2" in read or "READ2" in read.upper():
			self.right_reads.append(read)
		else:
			print("Warning: Paired read did not contain R1 or R2")
		

		
class UnpairedReadSample(Sample):

	def __init__(self, name):
		super().__init__(name)
		self.type = "unpaired"
			
			
		
	def getHisatString(self):
		hs = "-U "
		for i in range(len(self.reads)):
			read = self.reads[i]
			hs += str(read)
			if i < len(self.reads)-1:
				hs += ","	
		return hs
	
	
	
		
	def addRead(self, read):
		self.reads.append(read)
		


def extractGZIP(file):
	basename = file.split(".gz")[0]
	
	with gzip.open(file, "rt") as gzipped:
		with open(basename, "w") as outfile:
			for line in gzipped:
				outfile.write(line)
			
	


def canRun(exe):
	'''Runs "which" followed by the command and returns True if the returncode equals 0'''
	import subprocess
	try:
		subprocess.run(["which", exe], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		return True
	except:
		return False




def assureFolderEndsWithSlash(folder):
	if not folder.endswith("/") and not folder.endswith("\\"):
		folder+="/"
	return folder
	


def lookForGenome(folder="genome"):
	fileendings = [".fa", ".fasta", ".fas", ".fna"]
	
	folder = assureFolderEndsWithSlash(folder)
	
	for file in os.listdir(folder):
		for ending in fileendings:
			if file.lower().endswith(ending):
				return folder + file
	
	#could not find a fasta file. Maybe it is gzipped and needs to be extracted
	fileendings = [".fa.gz", ".fasta.gz", ".fna.gz"]
	for file in os.listdir(folder):
		for ending in fileendings:
			if file.lower().endswith(ending):
				print("Found a gzipped genome file. Extracting ...")
				extractGZIP(folder+file)
				return lookForGenome(folder)
	
	return None




def buildIndex(genome):
	basename = genome.split(".f")[0]
		
	cmd = "hisat2-build "
	cmd += str(genome)
	cmd += " "
	cmd += str(basename)
	print(cmd)
	os.system(cmd)
	

def mapToGenome(samples, genome, mappingFolder):
	genome_basename = genome.split(".f")[0]
	
	mappingFolder = assureFolderEndsWithSlash(mappingFolder)
	
	for sample in samples:
		cmd = "hisat2 -x " + str(genome_basename)
		cmd += " "
		cmd += "--threads "+threads + " "
		cmd += sample.getHisatString()
		cmd += " -S " + mappingFolder + sample.getName()
		cmd += ".sam"
		
		cmd += " --new-summary"
		cmd += " --summary-file "
		cmd += mappingFolder + sample.getName()+"_summary.txt"
		
		print(cmd)
		os.system(cmd)
		#immediately convert the SAM to BAM
		convertToBAM(mappingFolder)


def getGFFFile(folder="folder"):
	for file in os.listdir(folder):
		if file.lower().endswith("gff") or file.lower().endswith("gtf") or file.lower().endswith("gff3"):
			folder = assureFolderEndsWithSlash(folder)
			return folder+file
	
	return None



def runStringtie(gff_file, mapping_folder="mapping", outfolder="abundance/"):
	if not os.path.exists(outfolder):
		os.mkdir(outfolder)
		
	mapping_folder = assureFolderEndsWithSlash(mapping_folder)
	outfolder = assureFolderEndsWithSlash(outfolder)
	
	
	for file in os.listdir(mapping_folder):
		if not file.lower().endswith("bam") and not file.lower().endswith("sam"):
			continue
		name = file.split(".bam")[0].split(".sam")[0]
		
		cmd = "stringtie -e -B "
		cmd += "-p "+ threads + " "
		cmd+= mapping_folder+file
		cmd+= " -G "+gff_file
		cmd+= " -A "+outfolder + name+"/"+name + "_gene_expression.tsv"
		cmd+= " -o "+outfolder+name+"/"+name+".gtf"
		
		print(cmd)
		os.system(cmd)
	



def genomeIsAlreadyIndexed(folder="genome"):
	for file in os.listdir(folder):
		if file.endswith("ht2"):
			return True
	return False


def getAllReadsInFolder(folder="reads"):
	reads = []
	for file in os.listdir(folder):
		if isFastqFile(file):
			reads.append(file)
	return reads		
	
	

def getBasenameOfRead(read):
	basename = read.split(".fa")[0]
	basename = basename.split(".fq")[0]
	
	basename = basename.split("_R1")[0]
	basename = basename.split("_R2")[0]
	basename = basename.split("_Read1")[0]
	basename = basename.split("_Read2")[0]
	basename = basename.split("_read1")[0]
	basename = basename.split("_read2")[0]
	
	basename = basename.split("_L0")[0]
	return basename




def readAllMappingSummaries(folder="./mapping"):
	'''Reads all *_summary.txt files in a folder and returns a list
	 of the sample names, and a list of the text in the summary files.'''
	folder = assureFolderEndsWithSlash(folder)
	names = []
	summaries = []
	
	for file in os.listdir(folder):
		if not file.endswith("summary.txt"):
			continue
		
		with open(folder+file, "r") as summary:
			names.append(file.split("_summary")[0])
			text = ""
			for line in summary:
				text+=line
			summaries.append(text)

	return names,summaries






def extractNumberFromSummary(some_text_before_colon, text, extractPercentage=False):
	text = text[text.find(some_text_before_colon) : ]
	if extractPercentage:
		hit = re.search("\d+(\.?\d+?%)", text) #extracts all percentages
		hit = hit[0].replace("%", "")
		return float(hit)
	else:
		text = text.split(some_text_before_colon)[1]
		text = text.split(" (")[0]
		text = text.split("\n")[0]
		hit = re.search("\d+", text) #extracts all numbers
		text = int(hit[0])
	return text
	


def summarizePairedMapping(folder="./mapping"):
	try:
		import pandas as pd
	except:
		print("Warning: Could not create a summary of the mapping. Is pandas installed?")
		return
	
	folder = assureFolderEndsWithSlash(folder)
	
	names, summaries = readAllMappingSummaries(folder)
	data = dict()
	data["Sample"] = []
	data["Total pairs"] = []
	data["Aligned concordantly or discordantly 0 time"] = []
	data["Aligned concordantly or discordantly 0 time %"] = []
	data["Aligned concordantly 1 time"] = []
	data["Aligned concordantly 1 time %"] = []
	data["Aligned concordantly >1 times"] = []
	data["Aligned concordantly >1 times %"] = []
	data["Aligned discordantly 1 time"] = []
	data["Aligned discordantly 1 time %"] = []
	data["Total unpaired reads"] = []
	data["Aligned 0 time"] = []
	data["Aligned 0 time %"] = []
	data["Aligned 1 time"] = []
	data["Aligned 1 time %"] = []
	data["Aligned >1 times"] = []
	data["Aligned >1 times %"] = []
	data["Overall alignment rate"] = []

	for name, summary in zip(names, summaries):
		data["Sample"].append(name)
		data["Total pairs"].append(extractNumberFromSummary("Total pairs", summary))
		data["Aligned concordantly or discordantly 0 time"].append(extractNumberFromSummary("Aligned concordantly or discordantly 0 time", summary))
		data["Aligned concordantly or discordantly 0 time %"].append(extractNumberFromSummary("Aligned concordantly or discordantly 0 time", summary, True))
		data["Aligned concordantly 1 time"].append(extractNumberFromSummary("Aligned concordantly 1 time", summary))
		data["Aligned concordantly 1 time %"].append(extractNumberFromSummary("Aligned concordantly 1 time", summary, True))
		data["Aligned concordantly >1 times"].append(extractNumberFromSummary("Aligned concordantly >1 times", summary))
		data["Aligned concordantly >1 times %"].append(extractNumberFromSummary("Aligned concordantly >1 times", summary, True))
		data["Aligned discordantly 1 time"].append(extractNumberFromSummary("Aligned discordantly 1 time", summary))
		data["Aligned discordantly 1 time %"].append(extractNumberFromSummary("Aligned discordantly 1 time", summary, True))
		data["Total unpaired reads"].append(extractNumberFromSummary("Total unpaired reads", summary))
		data["Aligned 0 time"].append(extractNumberFromSummary("Aligned 0 time", summary))
		data["Aligned 0 time %"].append(extractNumberFromSummary("Aligned 0 time", summary, True))
		data["Aligned 1 time"].append(extractNumberFromSummary("Aligned 1 time", summary))
		data["Aligned 1 time %"].append(extractNumberFromSummary("Aligned 1 time", summary, True))
		data["Aligned >1 times"].append(extractNumberFromSummary("Aligned >1 times", summary))
		data["Aligned >1 times %"].append(extractNumberFromSummary("Aligned >1 times", summary, True))
		data["Overall alignment rate"].append(extractNumberFromSummary("Overall alignment rate", summary, True))
	
	df = pd.DataFrame(data)
	df.set_index("Sample", inplace=True)
	df.to_csv(folder+"mapping_summary.tsv",sep="\t")
	try:
		df.to_excel(folder+"mapping_summary.xlsx")
	except:
		pass
	
	return


def summarizeUnpairedMapping(folder="./mapping"):
	try:
		import pandas as pd
	except:
		print("Warning: Could not create a summary of the mapping. Is pandas installed?")
		return
	
	folder = assureFolderEndsWithSlash(folder)
	
	names, summaries = readAllMappingSummaries(folder)
	data = dict()
	data["Sample"] = []
	data["Total reads"] = []
	data["Aligned 0 time"] = []
	data["Aligned 0 time %"] = []
	data["Aligned 1 time"] = []
	data["Aligned 1 time %"] = []
	data["Aligned >1 times"] = []
	data["Aligned >1 times %"] = []
	data["Overall alignment rate"] = []


	for name, summary in zip(names, summaries):
		data["Sample"].append(name)
		data["Total reads"].append(extractNumberFromSummary("Total reads", summary))
		data["Aligned 0 time"].append(extractNumberFromSummary("Aligned 0 time", summary))
		data["Aligned 0 time %"].append(extractNumberFromSummary("Aligned 0 time", summary, True))
		data["Aligned 1 time"].append(extractNumberFromSummary("Aligned 1 time", summary))
		data["Aligned 1 time %"].append(extractNumberFromSummary("Aligned 1 time", summary, True))
		data["Aligned >1 times"].append(extractNumberFromSummary("Aligned >1 times", summary))
		data["Aligned >1 times %"].append(extractNumberFromSummary("Aligned >1 times", summary, True))
		data["Overall alignment rate"].append(extractNumberFromSummary("Overall alignment rate", summary, True))
	
	df = pd.DataFrame(data)
	df.set_index("Sample", inplace=True)
	df.to_csv(folder+"mapping_summary.tsv",sep="\t")
	try:
		df.to_excel(folder+"mapping_summary.xlsx")
	except:
		pass
	
	return
	
	
	

def mergeAbundances(abundance_folder):
	
	abundance_folder = assureFolderEndsWithSlash(abundance_folder)
	try:
		import pandas as pd
	except:
		print("Warning: Pandas is not installed. Cannot combine the FPKM's into a single file.")

	dataframe = None

	for subfolder in os.listdir(abundance_folder):
		if not os.path.isdir(abundance_folder+'/'+subfolder):
			continue

		df = pd.read_table(abundance_folder+subfolder+"/"+subfolder+"_gene_expression.tsv", index_col=0)
		df = df.drop(["TPM", "Coverage"], axis=1)

		if dataframe is None:
			dataframe = df
		else:
			dataframe = dataframe.join(df[["FPKM"]], how='outer')

		#lets rename the FPKM into the sample name
		renamed_columns = list(dataframe.columns)[0:-1]
		renamed_columns.append(subfolder)
		dataframe.columns = renamed_columns
	
	dataframe.to_csv(abundance_folder+"merged_FPKM.tsv", sep="\t")
	try:
		dataframe.to_excel(abundance_folder+"merged_FPKM.xlsx")
	except:
		pass
	
	return dataframe





def convertToBAM(folder):
	'''Converts the SAM files generated by hisat2 into BAM files using samtools'''
	folder = assureFolderEndsWithSlash(folder)
	
	for file in os.listdir(folder):
		if not file.lower().endswith("sam"):
			continue
		
		name = file.split(".sam")[0]
		
		cmd = "samtools sort -@ "+ str(threads)
		cmd += " -o "+folder+name+".bam " 
		cmd += folder+file
		
		print(cmd)
		returnCode = os.system(cmd)
		
		#remove the sam file if command ran without errors
		if returnCode == 0: 
			os.remove(folder+file)
		
		
	

def prepareReads(folder=".", arePaired=False):
	reads = getAllReadsInFolder(folder)
	sample_dict = dict()
	folder = assureFolderEndsWithSlash(folder)
	
	for r in reads:
		basename = getBasenameOfRead(r)
		#print(basename)
		entry = sample_dict.get(basename)
		
		if entry is None:
			entry = [folder+r]
		else:
			entry.append(folder+r)
		
		sample_dict[basename] = entry
	
	#print(sample_dict)
	
	samples = []
	for basename, reads in sample_dict.items():
		sample = None
		if arePaired:
			sample =  PairedReadSample(basename)
		else:
			sample =  UnpairedReadSample(basename)
		
		for read in reads:
			sample.addRead(read)
			
		samples.append(sample)
	
	return samples
	
	
def isFastqFile(file):
	if file.upper().endswith("FQ") or file.upper().endswith("FASTQ") or file.upper().endswith("FQ.GZ") or file.upper().endswith("FASTQ.GZ") :
		return True
	else:
		return False
		
		
		
def arePairedReads(folder="."):
	foundLeftRead = False
	foundRightRead = False
	
	for file in os.listdir(folder):
	
		if isFastqFile(file):
			if "R1" in file or "Read1" in file or "read1" in file:
				foundLeftRead = True
			elif "R2" in file or "Read2" in file or "read2" in file:
				foundRightRead = True

		if foundLeftRead and foundRightRead:
			return True
	return False



def checkIfAllPrerequisitesInstalled():
	print("Checking if the required software is installed...")
	found_hisat = canRun("hisat2")
	found_stringtie = canRun("stringtie")
	found_samtools = canRun("samtools")
	
	if not found_hisat:
		print("HISAT2 does not seem to be installed...")
		
	if not found_stringtie:
		print("stringtie does not seem to be installed...")
		
	if not found_samtools:
		print("samtools does not seem to be installed...")
	
	try:
		import pandas as pd
	except:
		print("Warning: Pandas is not installed. Your results won't be summarized.")


	
	if (not found_hisat) or (not found_samtools) or (not found_stringtie):
		print("")
		print("Some required software was not found.")
		print("Make sure it is installed, and that you're")
		print("running this script in the same environment.")
		print("Quitting...")
		quit()
	else:
		print("... no problems found.")
		print("")
		


print("HISAT2-pipline - Version 1.0.0 (2022/03) ")
print("")


checkIfAllPrerequisitesInstalled()


parser = argparse.ArgumentParser()
parser.add_argument("--skip_mapping", help="Skip mapping to genome", action="store_true")
parser.add_argument("--folder", help="Use this folder to look for the reads and genome folder", default=".")

args = parser.parse_args()

folder = args.folder
folder = assureFolderEndsWithSlash(folder)

threads = str(os.cpu_count())





PAIRED = arePairedReads(folder+"reads")
if PAIRED:
	paired_correct_answer = input("I found PAIRED reads in the folder. Is this correct? (yes/no) ")
	if "N" in paired_correct_answer.upper():
		PAIRED = not PAIRED
else:
	paired_correct_answer = input("I found UNPAIRED reads in the folder. Is this correct? (yes/no) ")
	if "N" in paired_correct_answer.upper():
		print("Please rename the file so that they contain _R1 and _R2 respectively (or _Read1 and _Read2). Otherwise I won't able to distinguish them.")
		print("Exiting")
		quit()	
	
if PAIRED:
	print("Running pipeline with PAIRED reads")
else:
	print("Running pipeline with UNPAIRED reads")

print()

samples = prepareReads(folder+"reads", PAIRED) 
print("Found the following samples:")
for sample in samples:
	print(sample)

print()

samples_correct = input("Is this correct? (yes/no) ")
if "n" in samples_correct.lower():
		print("Exiting")
		quit()	


print()
genome = lookForGenome(folder+"genome")
if genome is not None:
	if not genomeIsAlreadyIndexed(folder+"genome"):
		print("Building the genome index.")
		print()
		
		buildIndex(genome)
	else:
		reindex = input("I found a genome index in the genome folder. Do you want to skip building the index? (yes/no) ")
		if "n" in reindex.lower():
			print("Building the genome index.")
			print()
			buildIndex(genome)

else:
	print("Could not find a genome file in the genomes folder. Exiting")
	quit()
	
	
print()
if not os.path.exists(folder+"mapping"):
	os.mkdir(folder+"mapping")

if not args.skip_mapping:
	mapToGenome(samples, genome, folder+"mapping")
	if PAIRED:
		summarizePairedMapping(folder+"mapping")
	else:
		summarizeUnpairedMapping(folder+"mapping")
	
	
print("Converting SAM files to BAM files")	
convertToBAM(folder+"mapping")
print()


print("Estimating expression abundance with stringtie")
gff_file = getGFFFile(folder+"genome")
if gff_file is None:
	print("Could not find the GFF file for the genome. Are you sure it's in the genome folder?")
	print("Quitting..")
	quit()

runStringtie(gff_file, mapping_folder=folder+"mapping", outfolder=folder+"abundance")	
print()

mergeAbundances(folder+"abundance")
	
	
	
	
	
	
