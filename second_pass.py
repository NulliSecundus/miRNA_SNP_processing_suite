import click
import subprocess
import secrets
import time
import math
import io
from multiprocessing import Pool

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
rnaFileList = []
snpFileList = []
mirandaList = [] # Ordered parameters to iterate miranda
outputFileList = [] # List to hold names of temp miranda output files
noEnergy = False # miranda no energy option 
v = False
sigID = None
dir = None 
outputFileName = None 

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.option('--noenergy', is_flag=True, help='Flag for miranda -noenergy option')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(mirandafile, procsnpfile, mirnafile, output, noenergy, verbose):

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
	noEnergy = noenergy
	
	try:
		genSig(mirandafile, procsnpfile, mirnafile)
		loadsnp(procsnpfile)
		loadrna(mirnafile)
		'''
		print(len(procSnpArray))
		print(len(procSnpArray[0]))
		print(len(procSnpArray[1]))
		print(len(procSnpArray[143]))
		print(procSnpArray[0][1])
		print(procRnaArray[0])
		'''
		loadTopList(mirandafile)
		buildBottomList()
		iterateMiranda()
		print("Success")
	except Exception as error:
		toPrint = "Error: " + str(error)
		print(toPrint)
		'''
		print("Error")
		print(error)
		'''
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
	
	# Determine the number of entries in snp file 
	toRun = [
		"wc", 
		"-l",
		procSnpFasta
	]
	completedProcess = subprocess.run(toRun, stdout=subprocess.PIPE, encoding="utf-8")
	stdOutText = completedProcess.stdout
	textArray = stdOutText.split(" ")
	numLines = float(textArray[0])
	#print(numLines)
	snpSplit = int(math.ceil(math.sqrt(numLines/3)))
	
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
						
						if count%snpSplit == 0:
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
		
		printSnp()
	
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
				
	printRna()

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
				'''
				if count%100==0:
					break
				'''
				
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
	
	if v:
		toPrint = "Total data lines processed: " + str(count)
		print(toPrint)
	'''
	for entry in topList[0:10]:
		print(entry)
	
	for entry in topList[10000:10010]:
		print(entry)
	
	for entry in topList[-10:-1]:
		print(entry)
	'''
	validateTopListOrder()
	
	toPrint = "Length of topList: " + str(topListLen)
	#print(toPrint)
	
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
		bottomLine = [rsNum]
		
		for alleleList in entry[1:]:
			topAllele = alleleList[0]
			for checkAlleleList in snpLine[3:]:
				checkAllele = checkAlleleList[0]
				if checkAllele != topAllele:
					bottomAlleleList = [checkAllele] + alleleList[1:]
					bottomLine.append(bottomAlleleList)
		bottomList.append(bottomLine)
		'''
		print(bottomList)
		return
		'''
		
	if v:
		toPrint = "Length of bottomList: " + str(len(bottomList))
		print(toPrint)
		toPrint = bottomList[0]
		print(toPrint)

# Iteratively runs miranda on the list of SNP-miRNA pairs to be processed 
def iterateMiranda():
	global procSnpArray
	global procRnaArray
	
	# Clear memory of unused variables
	procSnpArray = None
	procRnaArray = None
	
	print("Processing list complete, running miranda")
	#return
	
	# Setup lists of files to iterate miranda over 
	for entry in snpFileList:
		tempMiranda = []
		for line in rnaFileList:
			tempMiranda.append([entry, line])
		mirandaList.append(tempMiranda)
	
	# Iteratively run miranda
	num = 1
	for entry in mirandaList:
		toPrint = "Running miranda on SNP part " + str(num)
		print(toPrint)
		with Pool() as p:
			p.map(runMiranda, entry)
		parseMiranda(num)
		appendOutput(outputFileName)
		num += 1
		
	# Delete all temp files and folder
	toPrint = "Clearing temporary files and folder"
	print(toPrint)
	toRun = ["rm", "-rf", dir.replace('/', '')]
	subprocess.run(toRun, check=True)
	
	# Add ending timestamp to output file 
	with open(outputFileName, 'a') as o:
		localtime = "# End: " + time.asctime(time.localtime(time.time()))
		print("{}".format(localtime), file=o)
	
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

# Utility function for handling each iteration of miranda 
def runMiranda(x):
	global dir
	global sigID
	
	tempSnpFile = x[0]
	tempRnaFile = x[1]
	
	# Get the file number for output naming
	textArray = tempRnaFile.replace('.fasta', '').split("_")
	fileNum = textArray[3]
	
	# Assign name of temp output file 
	tempOut = dir + sigID + "_out_" + fileNum + ".txt"
	
	# Run miranda
	if noEnergy:
		toRun = [
			"miranda", 
			tempRnaFile, 
			tempSnpFile, 
			"-sc",
			"0.0",
			"-noenergy",
			"-quiet",
			"-keyval",
			"-out",
			tempOut
		]
	else:
		toRun = [
			"miranda", 
			tempRnaFile, 
			tempSnpFile, 
			"-sc",
			"0.0",
			"-quiet",
			"-keyval",
			"-out",
			tempOut
		]
	subprocess.run(toRun, check=True)
	
# Utility function for parsing temp miranda outputs into final format 	
def parseMiranda(n):
	toPrint = "Parsing part " + str(n) + " miranda outputs"
	print(toPrint)
	
	for entry in outputFileList:
		output_container = [] # container for lines of final output
		container = [] # temporary container for comparison within a group
		
		with open(entry) as f:
			
			count = 0
			sumcount = 0
			
			for line in f:
				# for each line in the miranda output file
				if line[0:2]=='>h':
					# if line is a data line
					container.append(line.rstrip())
							
				elif line[0:2]=='>>':
					# if line is a summary line
					# end case for a grouping
					# copy top score from from container
					sumcount += 1
					topLine = ""
					
					if len(container)==1 :
						# if the group only has one data line 
						topLine = str(container[0])
					else:
						# if multiple data lines in grouping
						# determine top score line
						# send top score line to output_container
						topScore = -1
						
						for dataline in container:
							# for each data line in the grouping
							splitdline = dataline.split("\t")
							compare = float(splitdline[2])
							if compare > topScore :
								topScore = compare
								topLine = dataline
					
					# append topLine to output_container
					output_container.append(topLine)
						
					# reset the temporary container
					container = [] 
		
		cmpOut = entry.replace('.txt', '_cmp.txt')
		with open(cmpOut, 'w') as o:	
			for line in output_container:
				# write each topLine to the output file 
				print('{}'.format(line), file=o)

		# Delete temp miranda output file
		toRun = ["rm", entry]
		subprocess.run(toRun, check=True)
		
# Utility function for appending each compressed temp output to final
def appendOutput(out):
	for entry in outputFileList:
		cmpOut = entry.replace('.txt', '_cmp.txt')
		with open(cmpOut) as f:
			with open(out, 'a') as o:
				for line in f:
					print('{}'.format(line.replace('\n', '')), file=o)
		
		# Delete compressed temp output files
		toRun = ["rm", cmpOut]
		subprocess.run(toRun, check=True)
	
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
	
def searchTopSublist(rsNum, allele, rna, topStart, topStop):
	# Search in order
	index = 0
	for entry in topList[topStart:topStop]:
		cmpNum = entry[0]
		if rsNum < cmpNum:
			temp = [rsNum, [allele, rna]]
			topList.insert(topStart + index, temp)
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
	topList.insert(topStart + index, temp)
	
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
		searchTopSublist(rsNum, allele, rna, 0, topListLen)
		return
		
	# First check if rsNum is less than first entry in topList
	cmpNum = topList[0][0]
	if rsNum < cmpNum:
		temp = [rsNum, [allele, rna]]
		topList.insert(0, temp)
		return
		
	# Check if rsNum is greater than last entry in topList
	cmpNum = topList[-1][0]
	if rsNum > cmpNum:
		temp = [rsNum, [allele, rna]]
		topList.append(temp)
		return
	
	# Define increment value as sqrt of topListLen
	topSqrt = math.sqrt(topListLen)
	topListInc = int(math.ceil(topSqrt))
	topListIncSize = int(math.ceil(topSqrt))
		
	for x in range(topListIncSize):
		topStop = topStart + topListInc
		
		# If sublist is the last section of topList 
		# Special case search 
		if topStop > topListLen-1:
			searchTopSublist(rsNum, allele, rna, topStart, topListLen)
			return
			
		# Determine sublist boundaries
		topStartRs = topList[topStart][0]
		topStopRs = topList[topStop][0]
		
		# If rsNum falls between the boundaries
		if rsNum >= topStartRs and rsNum < topStopRs:
			# Search the sublist 
			#sublist = topList[topStart:topStop]
			searchTopSublist(rsNum, allele, rna, topStart, topStop)
			return
		
		# Otherwise, move to the next sublist 
		topStart = topStop
		
	'''
	# If rsNum is higher than all entries in list
	# Then append to end of topList 
	temp = [rsNum, [allele, rna]]
	topList.append(temp)
	'''
	
	print("Failed")
	print(rsNum)
	print(allele)
	print(rna)
	
# DEBUG: utility method for verify that topList is correctly sorted 
def validateTopListOrder():
	length = len(topList)
	
	for x in range(1, length-1):
		pre = topList[x-1][0]
		cur = topList[x][0]
		post = topList[x+1][0]
		if cur <= pre:
			print("Failed")
			toPrint = "Pre: " + str(pre)
			toPrint += " Cur: " + str(cur)
			toPrint += " Post: " + str(post)
			toPrint += " X: " + str(x)
			print(toPrint)
			return
		elif cur >= post:
			print("Failed")
			toPrint = "Pre: " + str(pre)
			toPrint += " Cur: " + str(cur)
			toPrint += " Post: " + str(post)
			toPrint += " X: " + str(x)
			print(toPrint)
			return
			
	print("Validated")
	
def printSnp():
	print("Printing SNP File")
	
	length = len(procSnpArray)
	
	# Loop through procSnpArray
	for n in range(length):
		sublist = procSnpArray[n]
		printSubSnp(sublist, n)
	
def printSubSnp(sublist, n):
	tempSnpFileName = dir + sigID + "_snp_" + str(n) + ".fasta"
	with open(tempSnpFileName, "w") as text_file:
		for entry in sublist[1:]:
			header = entry[1]
			alleleList = entry[3:]
			for allele in alleleList:
				subheader = header + allele[0]
				sequence = allele[1] + "\n"
				
				# Print to file 
				print("{}".format(subheader), file=text_file)
				print("{}".format(sequence), file=text_file)
	
	snpFileList.append(tempSnpFileName)
	
def printRna():
	print("Printing RNA File")
	
	rnaNum = len(procRnaArray)
	rnaSplit = int(math.ceil(rnaNum / 30))
	
	count = 0
	n = 0
	
	rnaSublist = []
	for entry in procRnaArray:
		rnaSublist.append(entry)
		count += 1
		if count%rnaSplit==0:
			printSubRna(rnaSublist, n)
			rnaSublist = []
			n += 1
			
	if count%rnaSplit!=0:
		printSubRna(rnaSublist, n)
		rnaSublist = []
		
def printSubRna(sublist, n):
	tempRnaFileName = dir + sigID + "_mirna_" + str(n) + ".fasta"
	outName = dir + sigID + "_out_" + str(n) + ".txt"
	with open(tempRnaFileName, "w") as text_file:
		for entry in sublist:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	rnaFileList.append(tempRnaFileName)
	outputFileList.append(outName)