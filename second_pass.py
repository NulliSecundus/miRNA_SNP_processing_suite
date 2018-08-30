import click
import subprocess
import secrets
import time
import math
import io
from multiprocessing import Pool
from operator import itemgetter

procSnpArray = []
procRnaArray = []
topList = []
topListLen = 0
bottomList = []
bottomRnaList = []
bulkRnaList = []
bottomSnpList = []
bulkSnpList = []
tempFileList = []
v = False
sigID = None
dir = None 
outputFileName = None 

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(mirandafile, procsnpfile, mirnafile, output, verbose):

	"""
	\b 
	Arguments:
		MIRANDAFILE - INPUT processed output from procfirstpass script 
		PROCSNPFILE - INPUT processed SNP file from procchrom script 
		MIRNAFILE - INPUT fasta-formatted mature miRNA file 
		OUTPUT - processed output from second pass of miRanda 
	
	Wrapper interface for running a miranda second pass."""

	global v
	global outputFileName
	
	v = verbose
	outputFileName = output
	
	try:
		genSig(mirandafile, procsnpfile, mirnafile)
		loadsnp(procsnpfile)
		loadrna(mirnafile)
		loadTopList(mirandafile)
		buildBottomList()
		"""
		processBottomList()
		iterateMiranda()
		"""
		print("Success")
	except:
		print("Error")
		return
	
# Generates a unique signature ID to be used in file naming
def genSig(mirandaFile, snpFile, mirnaFile):
	global sigID
	global dir 
	
	sigID = str(secrets.randbelow(999999999999))
	dir = "temp_" + sigID
	
	toRun = ["mkdir", dir]
	subprocess.run(toRun, check=True)
	
	dir += "/"
	
	# Create output file and add timestamp
	with open(outputFileName, 'w') as o:
		localtime = "# Start: " + time.asctime(time.localtime(time.time()))
		print("{}".format(localtime), file=o)
	
	# Create README info file in the temp folder
	infoFile = dir + "README" + ".txt"
	with open(infoFile, "w") as text_file:
		header = "Second Pass\nInput Parameters: "
		
		# Print to file 
		print("{}".format(header), file=text_file)
		print("{}".format(mirandaFile), file=text_file)
		print("{}".format(snpFile), file=text_file)
		print("{}".format(mirnaFile), file=text_file)
		print("{}".format(outputFileName), file=text_file)
	
# Loads the SNP sequence file into memory 
def loadsnp(procSnpFasta):
	global procSnpArray
	
	print("Loading SNP fasta file (may take a few minutes)")
	
	snpSubArray = []
	snpLine = []
	snpAllele = []
	headerSplit = None
	count = 0
	rsNum = 0
	alleleNum = 0
	header = ""
	snpName = ""
	allele = ""
	sequence = ""
	rsSet = -1
	
	with open(procSnpFasta) as f:
		for line in f:
			if line[0]==">":
				
				header = line
				header = header.replace('\n', '')
				headerSplit = header.split("|")
				
				snpName = headerSplit[0] + "|" + headerSplit[1] + "|" 
				snpName += headerSplit[2] + "|" + headerSplit[3] + "|"
				
				rsNum = int(headerSplit[2].replace("rs", ""))
				alleleNum = int(headerSplit[3])
				allele = headerSplit[4]
				
				if rsSet == rsNum:
					# Add to current entry 
					snpAllele = [allele, "seq"]
				else:
					if count != 0:
						# Store current SNP line 
						snpSubArray.append(snpLine)
						
						if count%3850 == 0:
							headerLine = [rsNum]
							snpSubArray.insert(0, headerLine)
							procSnpArray.append(snpSubArray)
							
							# Reset SNP sub-list
							snpSubArray = []
					
					# Create new SNP line entry
					snpLine = [rsNum, snpName, alleleNum]
					snpAllele = [allele, "seq"]
					rsSet = rsNum
					count += 1
				
			elif line[0]=="\n": 
				# End of sequence
				# populate each temp line with 
				# [snpName, rsEnd, alleleNum, allele, sequence, allele, sequence, etc...]
				
				snpAllele[1] = sequence
				snpLine.append(snpAllele)
					
			elif line[0]=="#":
				#Do nothing, it's a comment line
				pass
			else:
				# Sequence line 
				sequence = line
				sequence = sequence.replace('\n', '')
				
		snpSubArray.append(snpLine)
		headerLine = [rsNum+1]
		snpSubArray.insert(0, headerLine)
		procSnpArray.append(snpSubArray)
	
# Loads the miRNA file into memory 
def loadrna(miRNA):
	global procRnaArray
	
	print("Loading miRNA file")
	
	count = 0
	header = ""
	sequence = ""
	
	with open(miRNA) as f:
		for line in f:
			if line[0]==">":
				header = line
				header = header.replace('\n', '')
				
			else:
				# Sequence line 
				sequence = line
				sequence = sequence.replace('\n', '')
				
				temp = [header, sequence]
				procRnaArray.append(temp)
				
				count = count+1

# Populates the list of top hit pairs
def loadTopList(mirandaFile):
	global topList
	
	print('Loading list of top hits from miranda file')
	
	# Determine the number of lines in mirandaFile
	toRun = [
		"wc", 
		"-l", 
		mirandaFile, 
	]
	completedProcess = subprocess.run(toRun, stdout=subprocess.PIPE, encoding="utf-8")
	stdOutText = completedProcess.stdout
	textArray = stdOutText.split(" ")
	numLines = int(textArray[0])
	topSplit = int( math.floor(float(numLines) / 30.0 ) )
	
	# For each line in processed miranda output
	# Populate the top list with the miRNA and SNP rsNum pair
	
	subTopList = []
	count = 0
	
	with open(mirandaFile) as f:
		for line in f:
			if line[0]==">":
				lineSplit = line.split("\t")
				
				mirnaName = lineSplit[0]
				
				refName = lineSplit[1]
				refSplit = refName.split("|")
				rsNum = int(refSplit[2].replace("rs", ""))
				
				allele = refSplit[4]
				
				# Insert information into topList 
				insertTopList(rsNum, allele, mirnaName)
				
				count +=1
				if count%100==0:
					return
				
				'''
				rsEntry = searchTopRs(rsNum, allele, mirnaName)
				
				if rsEntry == None:
					temp = [rsNum, allele, mirnaName]
					topList.append(temp)
				
				
				temp = [mirnaName, rsNum, allele]
				subTopList.append(temp)
				count += 1
				
				if count%topSplit==0:
					topList.append(subTopList)
					subTopList = []
				'''
				
	for entry in topList[1000:1020]:
		print(entry)
		
	toPrint = "Length of topList: " + str(topListLen)
	print(toPrint)
	
# Populates the list of pairs to process with miranda
def buildBottomList():
	global topList
	global bottomList
	
	if v:
		toPrint = "Length of topList: " + str(len(topList))
		print(toPrint)
	
	print("Building list of entries to process (may take a few minutes)")
	
	# For each SNP-miRNA entry in the top list
	# Search the processed SNP list for alternative SNP alleles ID
	# Add the alternative alleles to the bottom list for processing
	for entry in topList:
		rsNum = entry[0]
		snpLine = snpSearch(rsNum)
		print(snpLine)
		return

# Populates the list of pairs to process with miranda
def buildBottomListOLD():
	global topList
	global bottomList
	global bulkRnaList
	global bulkSnpList
	
	if v:
		toPrint = "Length of topList: " + str(len(topList))
		print(toPrint)
	
	print("Building list of entries to process (may take a few minutes)")
	
	# For each SNP-miRNA entry in the top list
	# Search the processed SNP list for alternative SNP alleles ID
	# Add the alternative alleles to the bottom list for processing
	
	# Setup 30 entries in bottomList
	for x in range(30):
		bottomList.append(x)
	
	if v:
		toPrint = str(len(bottomList)) + " sublists in bottomList for multiprocessing support"
		print(toPrint)
	
	with Pool() as p:
		sublists = p.map(buildSubBottomList, bottomList)
	
	bottomList = []
	listNum = 1
	for list in sublists:
		'''
		print(list[0])
		print(listNum)
		'''
		
		bottomPairs = list[0]
		bulkRnaList.extend(list[1])
		bulkSnpList.extend(list[2])
		
		bottomPairs.insert(0, listNum)
		bottomList.append(bottomPairs)
		listNum += 1
	
# Sub-function to handle parallel processing of each topList section
def buildSubBottomList(n):
	
	count = 0
	sublist = []
	bulkRnaSublist = []
	bulkSnpSublist = []
	
	for line in topList[n]:
		
		mirna = line[0]
		rsNum = int(line[1])
		allele = line[2]
		
		snpLine = snpSearch(rsNum)
		
		if snpLine == None:
			toPrint = "Error in buildSubBottomList " + str(n)
			print(toPrint)
			pass
		
		snpName = snpLine[1]
		alleleNum = snpLine[2]
		
		for x in range(alleleNum):
			checkAllele = snpLine[3+x]
			if checkAllele[0] != allele:
				snpAlleleName = snpName + checkAllele[0]
				#sublist.append([mirna, mirnaSeq(mirna), snpAlleleName, checkAllele[1]])
				sublist.append([mirna, snpAlleleName])
				bulkRnaSublist.append([mirna])
				bulkSnpSublist.append([snpAlleleName])
		count += 1
	
	if v:
		toPrint = "Finished buildSubBottomList " + str(n) + ", sublist length " + str(len(sublist))
		print(toPrint)
		
	return [sublist, bulkRnaSublist, bulkSnpSublist]
	
def processBottomList():
	global bulkRnaList
	print("processing bottomList")
	
	toPrint = "BulkRnaList length: " + str(len(bulkRnaList))
	print(toPrint)
	
	count=0
	for mirna in bulkRnaList:
		if mirna not in bottomRnaList:
			bottomRnaList.append(mirna)
		count+=1
		
		if count%50000==0:
			print(count)
			toPrint = "BottomRnaList length: " + str(len(bottomRnaList))
			print(toPrint)
	
	toPrint = "BulkSnpList length: " + str(len(bulkSnpList))
	print(toPrint)
	
	count=0
	for snp in bulkSnpList:
		if snp not in bottomSnpList:
			bottomSnpList.append(snp)
		count+=1
		
		if count%50000==0:
			print(count)
			toPrint = "BottomSnpList length: " + str(len(bottomSnpList))
			print(toPrint)
				
	toPrint = "miRNA number: " + str(len(bottomRnaList))
	print(toPrint)
	toPrint = "SNP number: " + str(len(bottomSnpList))
	print(toPrint)

# Iteratively runs miranda on the list of SNP-miRNA pairs to be processed 
def iterateMiranda():
	global procSnpArray
	global procRnaArray
	
	# Clear memory of unused variables
	procSnpArray = None
	procRnaArray = None
	
	print("Processing list complete, running miranda")
	return
	
	'''
	# Create temp input files for each sublist in bottomList
	num = 1
	for list in bottomList:
		genInput(list, num)
		num += 1
	
	if v:
		toPrint = "Input files: " + dir
		print(toPrint)
		
		for entry in tempFileList:
			print(entry)
	'''
			
	# For each sublist in bottomList, run miranda in parallel 
	with Pool() as p:
		p.map(runMiranda, bottomList)
	
# Searches the procSnpArray and returns the line for the given rsNum
def snpSearch(rs):
	for subsection in procSnpArray:
		if rs < subsection[0][0]:
			for line in subsection:
				#print(line[0])
				if rs == line[0]:
					return line
	toPrint = "Error: could not find SNP " + str(rs)
	print(toPrint)
	return None
	
# Returns the sequence associated with the given miRNA name	
def mirnaSeq(mirnaName):
	for line in procRnaArray:
		mirnaCmp = line[0]
		mirnaCmp = mirnaCmp.split(" ")
		mirnaCmp = mirnaCmp[0]
		if mirnaCmp == mirnaName:
			return str(line[1])

# Runs miranda on the given SNP-miRNA pair from the bottom list  			
def runMiranda(x):
	# Create temp input text files for miRNA and SNP fasta seqs
	# Run miranda on each pair using temp input files
	# Output to temp output text file
	# Delete temp input text files
	# Open output text file, store line, delete output text file
	
	n = x[0]
	list = x[1:]
	
	tempOutputFile = dir + sigID + "_out_" + str(n) + ".txt"
	mirandaOutput = None
	count = 1
	
	for line in list:
		# For each line in the sublist
		# Create temp files containing single miRNA/SNP entries 
		# Then run miranda on these temp files 
		# Capture the output and append to output part file 
		tempRnaFileName = dir + sigID + "_mirna_" + str(n) + ".fasta"
		with open(tempRnaFileName, "w") as text_file:
			header = line[0]
			sequence = line[1]
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
			
		tempSnpFileName = dir + sigID + "_snp_" + str(n) + ".fasta"
		with open(tempSnpFileName, "w") as text_file:
			header = line[2]
			sequence = line[3]
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
		
		# Run miranda
		toRun = [
			"miranda", 
			tempRnaFileName, 
			tempSnpFileName, 
			"-sc",
			"0.0",
			"-noenergy",
			"-quiet",
			"-keyval"
		]
		completedProcess = subprocess.run(toRun, stdout=subprocess.PIPE, encoding="utf-8")
		mirandaText = completedProcess.stdout
		mirandaTextArray = mirandaText.split("\n")
		
		# Parse miranda text and append top hit to output file 
		parseMiranda(mirandaTextArray, tempOutputFile)
		
		'''
		if count%1000==0:
			print(count)
			
		count += 1
		'''
		
		'''
		# Delete temp input files
		toRun = ["rm", tempRnaFileName, tempSnpFileName]
		subprocess.run(toRun, check=True)
		'''
		
		'''
		with open(tempOutputFile, "a") as o:
			print("{}".format(mirandaText), file=o)
			
		break
		'''
	
def parseMiranda(text, out):
	container = [] # temporary container for comparison within a group
	
	for line in text:
		# for each line in the miranda output text
		if line[0:2]=='>h':
			# if line is a data line
			container.append(line.rstrip())
					
		elif line[0:2]=='>>':
			# if line is a summary line
			# end case for a grouping
			# copy top score from from container
			topLine = ""
			
			if len(container)==1 :
				# if the group only has one data line 
				topLine = str(container[0])
			else:
				# if multiple data lines in grouping
				# determine top score line
				topScore = -1
				
				for dataline in container:
					# for each data line in the grouping
					splitdline = dataline.split("\t")
					compare = float(splitdline[2])
					if compare > topScore :
						topScore = compare
						topLine = dataline
				
			# reset the temporary container
			container = [] 
	
	with open(out, 'a') as o:	
		# write topLine to the output file 
		print('{}'.format(topLine), file=o)
	
def genInput(sublist, n):
	tempSnpFileName = dir + sigID + "_snp_" + str(n) + ".fasta"
	tempRnaFileName = dir + sigID + "_mirna_" + str(n) + ".fasta"
	
	with open(tempRnaFileName, "w") as r, open(tempSnpFileName, "w") as s:
		for entry in sublist:
			rnaHeader = entry[0]
			rnaSequence = entry[1] + "\n"
			snpHeader = entry[2]
			snpSequence = entry[3] + "\n"
			
			# Print to file 
			print("{}".format(rnaHeader), file=r)
			print("{}".format(rnaSequence), file=r)
			print("{}".format(snpHeader), file=s)
			print("{}".format(snpSequence), file=s)
	
	tempFileList.append([tempRnaFileName, tempSnpFileName])
	
def searchBottomRna(rna):
	# True means the miRNA was not found
	# False means the miRNA was found
	
	# if the list is empty return False
	if len(bottomRnaList)==0:
		return True
		
	# if the list is not empty then perform search
	for entry in bottomRnaList:
		if entry == rna:
			# return True if found
			return False
			
		# return False if not found 
	return True
	
def searchTopSublist(rsNum, allele, rna, sublist):
	# Search in order
	index = 0
	for entry in sublist:
		cmpNum = entry[0]
		if rsNum < cmpNum:
			temp = [rsNum, [allele, rna]]
			sublist.insert(index, temp)
			return
		elif rsNum == cmpNum:
			# If a matching rsNum is found
			# Check the allele(s)
			for alleleList in entry[1:]:
				cmpAllele = alleleList[0]
				if cmpAllele == allele:
					alleleList.append(rna)
					return
			temp = [allele, rna]
			entry.append(temp)
			return
		index += 1
		
	# If new rsNum is higher than all entries in list 
	# Then append to end of the list 
	temp = [rsNum, [allele, rna]]
	sublist.append(temp)
	
def insertTopList(rsNum, allele, rna):
	global topListLen
	topListLen = len(topList)
	
	topStart = 0
	topStop = 0
	
	if topListLen==0:
		temp = [rsNum, [allele, rna]]
		topList.append(temp)
		return
		
	elif topListLen < 10:
		searchTopSublist(rsNum, allele, rna, topList)
		return
		
	# Define increment value as sqrt of topListLen
	topSqrt = math.sqrt(topListLen)
	topListInc = int(math.floor(topSqrt))
	topListIncSize = int(math.ceil(topSqrt))
		
	for x in range(topListIncSize):
		topStop = topStart + topListInc
		if topStop > (topListLen-1):
			topStop = topListLen-1
			
		# Determine sublist boundaries
		topStartRs = topList[topStart][0]
		topStopRs = topList[topStop][0]
		
		# If rsNum falls between the boundaries
		if rsNum >= topStartRs and rsNum <= topStopRs:
			# Search the sublist 
			sublist = topList[topStart:topStop]
			searchTopSublist(rsNum, allele, rna, sublist)
			return
			
		# If rsNum is higher than all entries in list
		elif rsNum > topStopRs:
			temp = [rsNum, [allele, rna]]
			topList.append(temp)
			return
		
		# Otherwise, move to the next sublist 
		topStart = topStop
		
	print("Failed")