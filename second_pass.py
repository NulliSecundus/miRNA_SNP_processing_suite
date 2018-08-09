import click
import subprocess
import secrets
import time
import math
from multiprocessing import Pool

procSnpArray = []
snpList = [] # List to hold SNP entries
snpFileList = [] # List to hold names of SNP temp files 
rnaList = [] # List to hold miRNA entries
rnaFileList = [] # List to hold names of miRNA temp files 
topList = []
bottomList = []
mirandaList = [] # Ordered parameters to iterate miranda
outputFileList = [] # List to hold names of temp miranda output files
outputFileName = None
tempRestrictFileName = None 
sigID = None # Unique signature ID for naming temp files and folder 
dir = None # Directory for temp file storage 
sc = None # Input score threshold (float) 
snpSplit = 2000000 # Number of SNP entries per subsection
rnaSplit = 80 # Number of miRNA entries per subsection
noEnergy = False # miranda no energy option 

@click.command()
@click.argument('mirandafile')
@click.argument('snpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.argument('score')
@click.option('-snp', default=2000000, help='Number of SNP entries per subsection\nDefault 2000000')
@click.option('-rna', default=80, help='Number of miRNA entries per subsection\nDefault 80')
@click.option('-stop', default=0, help='Limits run to the specified number of SNP subsections')
@click.option('--noenergy', is_flag=True, help='Flag for miranda -noenergy option')
def cli(mirandafile, snpfile, mirnafile, output, score, snp, rna, stop, noenergy):
	global sc 
	global snpSplit
	global rnaSplit
	global noEnergy
	
	sc = float(score)
	outputFileName = output
	snpSplit = snp 
	rnaSplit = rna 
	noEnergy = noenergy
	
	try:
		genSig(mirandafile, procsnpfile, mirnafile)
		printsnp(snpfile)
		printrna(mirnafile)
		loadsnp(snpfile)
		loadTopList(mirandafile)
		buildBottomList()
		iterateMiranda(stop)
		print("Success")
	except:
		print("Error")
		return
	
# Generates a unique signature ID to be used in file naming
def genSig(mirandaFile, snpFile, mirnaFile):
	global sigID
	global dir 
	global tempRestrictFileName
	
	sigID = str(secrets.randbelow(999999999999))
	dir = "temp_" + sigID
	
	toRun = ["mkdir", dir]
	subprocess.run(toRun, check=True)
	
	dir += "/"
	tempRestrictFileName = dir + sigID + "_restrict.txt"
	
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

# Loads the SNP file into memory and splits it every 2 million entries
# Outputs each split section as a temp file for miranda 	
def printsnp(snpFile):
	global snpList
	
	print("Loading SNP fasta file")
	
	snpBlock = None
	fileNum = 1
	count = 0
	
	with open(snpFile) as f:
		for line in f:
			if line[0]==">":
				snpBlock = [line.replace('\n', ''), "seq"]
			elif line[0]=="\n": 	
				if snpBlock != None:
					snpList.append(snpBlock)
					count += 1
					if count%snpSplit==0:
						outputSnp(fileNum)
						fileNum += 1
			else:
				snpBlock[1] = line.replace('\n', '')
	
	# Output any remaining SNP entries
	if count%snpSplit!=0:
		outputSnp(fileNum)
	
# Loads the miRNA file into memory and splits it every 80 entries 
# Outputs each split section as a temp file for miranda 
def printrna(mirnaFile):
	global rnaList
	global outputFileList
	
	print("Loading miRNA fasta file")
	
	rnaBlock = None
	fileNum = 1
	count = 0
	
	with open(mirnaFile) as f:
		for line in f:
			if line[0]==">":
				rnaBlock = [line.replace('\n', ''), "seq"]
			else:
				rnaBlock[1] = line.replace('\n', '')
				rnaList.append(rnaBlock)
				count += 1
				if count%rnaSplit==0:
					outputRna(fileNum)
					fileNum += 1
				
	# Output any remaining miRNA entries
	if count%rnaSplit!=0:
		outputRna(fileNum)
		
	# Populate list of output files based on miRNA input files
	for n in range(fileNum):
		outName = dir + sigID + "_out_" + str(n+1) + ".txt"
		outputFileList.append(outName)
		
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
				
				temp = [mirnaName, rsNum, allele]
				subTopList.append(temp)
				count += 1
				
				if count%topSplit==0:
					topList.append(subTopList)
					subTopList = []
					
# Populates the list of pairs to process with miranda
def buildBottomList():
	global topList
	global bottomList
	
	if v:
		toPrint = "topList dimensions " + str(len(topList)) + " x " + str(len(topList[0]))
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
		#list.insert(0, listNum)
		
		#bottomList.append(list)
		#listNum += 1
		
		genRestrict(list)
			
	
# Sub-function to handle parallel processing of each topList section
def buildSubBottomList(n):
	count = 0
	sublist = []
	
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
				sublist.append([mirna, "seq", snpAlleleName, checkAllele[1]])
		count += 1
	
	if v:
		toPrint = "Finished buildSubBottomList " + str(n) + ", sublist length " + str(len(sublist))
		print(toPrint)
		
	return sublist
	
def iterateMiranda(stop):
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
		
		if (stop!=0) and (num%stop==0):
			break
		
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
	
# Utility function for printing current snpList to file
def outputSnp(n):
	global snpList
	global sigID
	global dir
	
	tempSnpFileName = dir + sigID + "_snp_" + str(n) + ".fasta"
	with open(tempSnpFileName, "w") as text_file:
		for entry in snpList:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	snpFileList.append(tempSnpFileName)
	snpList = []
	
# Utility function for printing current rnaList to file
def outputRna(n):
	global rnaList
	global sigID
	global dir
	
	tempRnaFileName = dir + sigID + "_mirna_" + str(n) + ".fasta"
	with open(tempRnaFileName, "w") as text_file:
		for entry in rnaList:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	rnaFileList.append(tempRnaFileName)
	rnaList = []
	
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
			str(sc),
			"-noenergy",
			"-quiet",
			"-keyval",
			"-restrict",
			tempRestrictFileName,
			"-out",
			tempOut
		]
	else:
		toRun = [
			"miranda", 
			tempRnaFile, 
			tempSnpFile, 
			"-sc",
			str(sc),
			"-quiet",
			"-keyval",
			"-restrict",
			tempRestrictFileName,
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
	
def genRestrict(sublist):
	with open(tempRestrictFileName, "a") as res:
		
		for entry in sublist:
			rnaHeader = entry[0]
			snpHeader = entry[2]
			restrict = rnaHeader.replace(">", "") + "\t" + snpHeader.replace(">", "")
			
			# Print to restrict file
			print("{}".format(restrict), file=res)