import click
import subprocess
import secrets

procSnpArray = []
procRnaArray = []
topList = []
bottomList = []

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(mirandafile, procsnpfile, mirnafile, output, verbose):
	try:
		loadsnp(procsnpfile)
		loadrna(mirnafile)
		loadTopList(mirandafile)
		buildBottomList()
		addSequences()
		iterateMiranda(output)
		print("Success")
	except:
		print("Error")
		return
	
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
						
						if count%5000 == 0:
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
				
				'''
				if count%1000000 == 0:
					print(count)
					
					for line in procSnpArray:
						print(line)
					
					return
					'''
					
			elif line[0]=="#":
				#Do nothing, it's a comment line
				pass
			else:
				# Sequence line 
				sequence = line
				sequence = sequence.replace('\n', '')
				
		snpSubArray.append(snpLine)
		headerLine = [rsNum]
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
				
				'''
				if count%100000 == 0:
					print(count)
					
					for line in procRnaArray:
						print(line)
					
					return
				'''
				
				count = count+1

# Populates the list of top hit pairs
def loadTopList(mirandaFile):
	global topList
	
	print('Loading list of top hits from miranda file')
	
	# For each line in processed miranda output
	# Populate the top list with the miRNA and SNP rsNum pair 
	
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
				topList.append(temp)

# Populates the list of pairs to process with miranda
def buildBottomList():
	global topList
	global bottomList
	
	print(len(topList))
	print("Building list of entries to process (may take a few minutes)")
	
	# For each SNP-miRNA entry in the top list
	# Search the processed SNP list for alternative SNP alleles ID
	# Add the alternative alleles to the bottom list for processing
	
	count = 0
	
	for line in topList:
		if count%100000==0:
			print(count)
		
		mirna = line[0]
		rsNum = int(line[1])
		allele = line[2]
		
		snpLine = snpSearch(rsNum)
		#print(snpLine)
		snpName = snpLine[1]
		alleleNum = snpLine[2]
		for x in range(alleleNum):
			checkAllele = snpLine[3+x]
			if checkAllele[0] != allele:
				snpAlleleName = snpName + checkAllele[0]
				bottomList.append([mirna, mirnaSeq(mirna), snpAlleleName, checkAllele[1]])
		count += 1
	
	print(bottomList[0])
	return

# Adds sequences to each identification label in the bottom list 
def addSequences():
	print('Loading sequences into bottom list (may take a few minutes)')
	mirnaName = ""
	snpName = ""
	count = 0
	
	for line in reprocessList:
		mirnaName = line[0]
		snpName = line[1]
		line.insert(1, mirnaSeq(mirnaName))
		line[2] = snpSeq(snpName)
		count += 1
		
		'''
		if count%100000==0:
			print(count)
		'''

# Iteratively runs miranda on the list of SNP-miRNA pairs to be processed 
def iterateMiranda(outputFile):
	global procSnpArray
	global procRnaArray
	
	# Clear memory of unused variables
	procSnpArray = None
	procRnaArray = None
	
	sig = str(secrets.randbelow(999999999999))
	count = 0
	
	'''
	print(len(reprocessList))
	'''
	
	print("Success: reprocess list complete, running miranda")
	
	# Iterate through list of SNP-miRNA pairs that need to be reprocessed 
	# Add score line to condensed final output file
	
	#outputFile = "repass1_sc206_chr1.txt"
	errorFile = outputFile.strip(".txt") + "_error_log.txt"
	with open(outputFile, "a") as final_output:
		with open(errorFile, "a") as error_log:
			for line in reprocessList:
				scoreLine = runMiranda(line, sig)
				
				# Print to file 
				if scoreLine != None:
					print("{}".format(scoreLine), file=final_output)
				else:
					print("{}".format(line), file=error_log)
				
				count += 1
				'''
				if count%100000==0:
					print(count)
				'''

# Returns the sequence associated with the given SNP name
def snpSeq(snpName):
	nameSplit = snpName.split("|")
	rsText = nameSplit[2]
	rsNum = int(rsText[2:])

	for line in procSnpArray:
		rsStart = line[0][0]
		rsEnd = line[0][1]

		if ((rsNum > rsStart) and (rsNum <= rsEnd)):

			for entry in line[1:]:

				cmpRsNum = entry[1]
				if cmpRsNum == rsNum:
					'''
					seqArray = [[label1, seq1]
								[label2, seq2]
								[label3, seq3]]
								...
								'''
					seqArray = []
					snpIndex = 1
					for n in range(entry[2]):
						label = ">" + snpName + "|" + str(entry[2])
						label += "|" + entry[snpIndex+2]
						temp = [label, entry[snpIndex+3]]
						
						seqArray.append(temp)
						snpIndex += 2
					return seqArray
					
	print("Failed to locate SNP Sequence")
	print(snpName)
	
# Searches the procSnpArray and returns the line for the given rsNum
def snpSearch(rs):
	for subsection in procSnpArray:
		if rs < subsection[0][0]:
			for line in subsection:
				#print(line[0])
				if rs == line[0]:
					return line
	print("Error: could not find SNP")
	
# Returns the sequence associated with the given miRNA name	
def mirnaSeq(mirnaName):
	for line in procRnaArray:
		mirnaCmp = line[0]
		mirnaCmp = mirnaCmp.split(" ")
		mirnaCmp = mirnaCmp[0]
		if mirnaCmp == mirnaName:
			return str(line[1])

# Runs miranda on the given SNP-miRNA pair from the bottom list  			
def runMiranda(reprocessLine, sig):
	# Create temp input text files for miRNA and SNP fasta seqs
	# Run miranda on each pair using temp input files
	# Output to temp output text file
	# Delete temp input text files
	# Open output text file, store line, delete output text file
	
	toReturn = None
		
	tempmirna = "temp_mirna_" + sig + "_input.fasta"
	with open(tempmirna, "w") as text_file:
		header = reprocessLine[0]
		sequence = reprocessLine[1]
		
		# Print to file 
		print("{}".format(header), file=text_file)
		print("{}".format(sequence), file=text_file)
		
	tempsnp = "temp_snp_" + sig + "_input.fasta"
	with open(tempsnp, "w") as text_file:
		snpArray = reprocessLine[2]
		for entry in snpArray:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	# Run miranda
	toRun = [
		"miranda", 
		tempmirna, 
		tempsnp, 
		"-sc",
		"206.0",
		"-noenergy",
		"-quiet",
		"-keyval"
	]
	completedProcess = subprocess.run(toRun, stdout=subprocess.PIPE, encoding="utf-8")
	mirandaText = completedProcess.stdout
	mirandaTextArray = mirandaText.split("\n")
	
	for line in mirandaTextArray:
		if line[0:2]=='>h':
			toReturn = line.rstrip()
	
	'''
	if toReturn==None:
		print(mirandaText)
	'''
	
	# Delete temp input files
	toRun = ["rm", tempmirna, tempsnp]
	subprocess.run(toRun, check=True)
	
	# Return the score line
	return toReturn