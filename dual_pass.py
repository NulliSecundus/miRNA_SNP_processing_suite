import click
import subprocess
import secrets
import time
import math
import io
from multiprocessing import Pool

procSnpArray = [] # List to hold SNP entries
procRnaArray = [] # List to hold miRNA entries
tempFileList = []
rnaFileList = [] # List to hold names of miRNA temp files 
snpFileList = [] # List to hold names of SNP temp files 
mirandaList = [] # Ordered parameters to iterate miranda
outputFileList = [] # List to hold names of temp miranda output files
noEnergy = False # miranda no energy option 
v = False
sigID = None # Unique signature ID for naming temp files and folder 
dir = None # Directory for temp file storage 
outputFileName = None
upperThreshold = 0.0 # Upper score threshold cutoff 
lowerThreshold = 0.0 # Lower score threshold cutoff

@click.command()
@click.argument('procsnpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.option('-ut', default=102.0, help='Upper score threshold cutoff\nDefault 102.0')
@click.option('-lt', default=72.0, help='Lower score threshold cutoff\nDefault 72.0')
@click.option('--noenergy', is_flag=True, help='Flag for miranda -noenergy option')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(procsnpfile, mirnafile, output, ut, lt, noenergy, verbose):

	"""
	\b 
	Arguments:
		PROCSNPFILE - INPUT processed SNP file from procchrom script 
		MIRNAFILE - INPUT fasta-formatted mature miRNA file 
		OUTPUT - processed output from dual pass of miranda
	
	Wrapper interface for running a miranda dual pass."""

	global v
	global outputFileName
	global upperThreshold
	global lowerThreshold
	
	v = verbose
	outputFileName = output
	noEnergy = noenergy
	upperThreshold = ut
	lowerThreshold = lt
	
	try:
		genSig(procsnpfile, mirnafile)
		loadsnp(procsnpfile)
		loadrna(mirnafile)
		iterateMiranda()
		print("Success")
	except Exception as error:
		toPrint = "Error: " + str(error)
		print(toPrint)
		return
	
# Generates a unique signature ID to be used in file naming
def genSig(snpFile, mirnaFile):
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
		thresholdText = "# Upper score threshold: " + str(upperThreshold)
		thresholdText += " Lower score threshold: " + str(lowerThreshold)
		print("{}".format(localtime), file=o)
		print("{}".format(thresholdText), file=o)
	
	# Create README info file in the temp folder
	infoFile = dir + "README" + ".txt"
	with open(infoFile, "w") as text_file:
		header = "Dual Pass\nInput Parameters: "
		ut = "Upper score threshold: " + str(upperThreshold)
		lt = "Lower score threshold: " + str(lowerThreshold)
		
		# Print to file 
		print("{}".format(header), file=text_file)
		print("{}".format(snpFile), file=text_file)
		print("{}".format(mirnaFile), file=text_file)
		print("{}".format(outputFileName), file=text_file)
		print("{}".format(ut), file=text_file)
		print("{}".format(lt), file=text_file)
	
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
	upperContainer = []
	lowerContainer = []
	cmpContainer = []
	cmpRs = 0
	cmpAlleleNum = 10
	
	for entry in outputFileList:
		cmpOut = entry.replace('.txt', '_cmp.txt')
		with open(cmpOut) as f, open(out, 'a') as o:
			for line in f:
				lineSplit = line.split("\t")
				
				mirnaName = lineSplit[0]
				refName = lineSplit[1]
				refScore = float(lineSplit[2])
				refEn = float(lineSplit[3])
				
				refSplit = refName.split("|")
				rsNum = int(refSplit[2].replace("rs", ""))
				alleleNum = int(refSplit[3])
				allele = refSplit[4]
				
				if rsNum == cmpRs :
					# add to container
					tempLine = [line, refScore]
					cmpContainer.append(tempLine)
				else:
					# if rsNum is diff than cmpRs
					# count num in container to verify all present
					if len(cmpContainer) == cmpAlleleNum:
						for entry in cmpContainer:
							# pass lines above upperThreshold to upperContainer
							if entry[1] >= upperThreshold:
								upperContainer.append(entry[0])
							# pass lines below lowerThreshold to lowerContainer
							if entry[1] <= lowerThreshold:
								lowerContainer.append(entry[0])
					if (len(upperContainer) > 0) and (len(lowerContainer) > 0):
						# if upperContainer and lowerContainer both contain entries		
						# write upperContainer and lowerContainer to output
						for upperLine in upperContainer:
							print('{}'.format(upperLine.replace('\n', '')), file=o)
						for lowerLine in lowerContainer:
							print('{}'.format(lowerLine.replace('\n', '')), file=o)
						
					# start new container 
					cmpContainer = []
					tempLine = [line, refScore]
					cmpContainer.append(tempLine)
					
					cmpRs = rsNum
					cmpAlleleNum = alleleNum
					upperContainer = []
					lowerContainer = []
		
		# Delete compressed temp output files
		toRun = ["rm", cmpOut]
		subprocess.run(toRun, check=True)
	
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