import click
import subprocess

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('mirna')
@click.argument('output')
def cli(mirandafile, procsnpfile, mirna, output):
	try:
		loadsnp(procsnpfile)
		loadrna(mirna)
		buildReprocList(mirandafile)
		addSequences()
		#return
		iterateMiranda()
	except:
		print("Error")
		return

if __name__ == '__main__':
	pass

snpInfo = []
mirnaInfo = []
reprocessList = []
	
def loadsnp(procSnpFasta):
	
	try:
		print("Loading SNP fasta file")
		
		snpSubInfo = []
		count = 0
		header = ""
		snpName = ""
		allele = ""
		alleleNum = 0
		alleleIndex = 0
		sequence = ""
		rsStart = 0
		rsEnd = 0
		rsSet = 0
		temp = []
		
		with open(procSnpFasta) as f:
			for line in f:
				if line[0]==">":
					
					if alleleIndex==alleleNum:
						alleleIndex = 0
					
					header = line
					header = header.replace('\n', '')
					
					snpName = header.split(" ")
					snpName = snpName[0][1:]
                    
					headerSplt = header.split("|")
					
					rsNum = headerSplt[2].split(" ")
					rsNum = rsNum[0][2:]
					rsEnd = int(rsNum)
					
					alleles = headerSplt[8]
					alleles = alleles[8:]
					alleles = alleles.replace('\"','')
					alleles = alleles.split('/')
					alleleNum = len(alleles)
					allele = alleles[alleleIndex]
					
					'''
					snpName = snpName + "|" + str(alleleNum)
					snpName = snpName + "|" + alleles[alleleIndex]
					'''
					
					alleleIndex = alleleIndex + 1
					
				elif line[0]=="\n": 
					# End of sequence
					# populate each temp line with 
					# [snpName, rsEnd, alleleNum, allele, sequence, allele, sequence, etc...]
					
					if rsSet == rsEnd:
						temp.extend([allele, sequence])
					else:
						if count != 0:
							snpSubInfo.append(temp)
							'''
							if rsSet==960618749: 
								print(temp
							'''
						temp = [snpName, rsEnd, alleleNum, allele, sequence]
						
						count = count+1
						if count%5000 == 0:
							headerLine = [rsStart, rsSet]
							snpSubInfo.insert(0, headerLine)
							snpInfo.append(snpSubInfo)
							
							# Reset SNP sub-list and rsStart number
							snpSubInfo = []
							rsStart = rsSet
						
						rsSet = rsEnd
					'''
					if count%1000000 == 0:
						print(count)
						
						for line in snpInfo:
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
					
			snpSubInfo.append(temp)
			headerLine = [rsStart, rsEnd]
			snpSubInfo.insert(0, headerLine)
			snpInfo.append(snpSubInfo)
	except:
		print('Could not parse processed snp fasta file')
	
def loadrna(miRNA):
	
	try:
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
					mirnaInfo.append(temp)
					
					'''
					if count%100000 == 0:
						print(count)
						
						for line in mirnaInfo:
							print(line)
						
						return
					'''
					
					count = count+1

	except:
		print('Could not parse miRNA file')

def buildReprocList(mirandaFile):
	try: 
		print('Building reprocess list from miranda file')
		
		# For each line in cleaned miranda output
		# Populate temp container with same-name entries
		# When new name appears, count num in temp container
		# Search snpInfo for the entry, compare to temp container num
		# If (available - output) > 0, then add SNP-miRNA pair to reprocess list
		
		#tempCont = []
		count = 0
		current = ""
		alleleCount = 1
		
		with open(mirandaFile) as f:
			for line in f:
				if line[0]==">":
					lineEdit = line.split("\t")
					mirnaName = lineEdit[0]
					refName = lineEdit[1]
					
					if refName != current:
						if (current != "") and (alleleCount==1):
							#checkAlleleCount(current, alleleCount, mirnaName)
							temp = [mirnaName, current]
							reprocessList.append(temp)
						alleleCount = 1
						current = refName
					else:
						alleleCount += 1
					
					count = count+1
				else:
					pass
				'''	
				if count%1000000 == 0:
					print(count)
					
					for entry in reprocessList:
						print(entry)
					
					return
					'''
				
	except:
		print('Could not build reprocess list from miranda file')

def addSequences():
	print('Loading sequences into reprocess list')
	mirnaName = ""
	snpName = ""
	count = 0
	
	for line in reprocessList:
		mirnaName = line[0]
		snpName = line[1]
		line.insert(1, mirnaSeq(mirnaName))
		line[2] = snpSeq(snpName)
		count += 1
		
		if count%100000==0:
			#if count%10==0:
			print(count)
			#return
			
	'''
	count=0
	outputFile = "chr1_restrict_file.txt"
	with open(outputFile, "a") as final_output:
		for line in reprocessList:	
			mirnaName = line[0]
			snpArray = line[2]
			
			# Print to file 
			for entry in snpArray:
				toPrint = mirnaName + "\t" + entry[0]
				print("{}".format(toPrint), file=final_output)
			
			count += 1
			
			if count%100000==0:
				print(count)
	'''
		
def iterateMiranda():
	try: 
		# Clear memory of unused variables
		snpInfo = None
		mirnaInfo = None
		
		count = 0
		
		'''
		print(reprocessList[0])
		print(len(reprocessList))
		'''
		print("Success: reprocess list complete, running miranda")
		
		# Iterate through list of SNP-miRNA pairs that need to be reprocessed 
		# Add score line to condensed final output file
		
		outputFile = "condensed_chr1_pass1.txt"
		with open(outputFile, "a") as final_output:
			for line in reprocessList:
				scoreLine = runMiranda(line)
				
				# Print to file 
				if scoreLine != None:
					print("{}".format(scoreLine), file=final_output)
				else:
					print("Error")
					print(line)
				
				count += 1
				
				if count%500==0:
					print(count)

	except:
		print("Failed to run miranda on reprocess list")
		
# DEPRECIATED: checks the available vs. found allele num for a given SNP seq 
def checkAlleleCount(name, num, mirna):
	#print(name + " " + str(num))
	#startSig = False 
	#temp = [mirna, mirnaSeq(mirna)]
	
	for line in snpInfo:
		strcmp = line[0]
		if strcmp == name:
			if (line[1] - num) > 0:
				'''
				if not(startSig):
					temp = [mirna, mirnaSeq(mirna)]

				startSig = True 
				'''
				#print("re-process " + name)
				#refSeqName = line[0] + "|" + str(line[1]) + "|" + line[2]
				# Add in alleles 
				#temp += [refSeqName, line[3]]
				temp = [mirna, strcmp]
				reprocessList.append(temp)
				return True
			else:
				return False
		'''
		else: 
			if startSig:
				reprocessList.append(temp)
				return
		'''
	return False

# Returns the sequence associated with the given SNP name
def snpSeq(snpName):
	nameSplit = snpName.split("|")
	rsText = nameSplit[2]
	rsNum = int(rsText[2:])
	#print(rsNum)

	for line in snpInfo:
		rsStart = line[0][0]
		rsEnd = line[0][1]

		if ((rsNum > rsStart) and (rsNum <= rsEnd)):

			for entry in line[1:]:

				cmpRsNum = entry[1]
				if cmpRsNum == rsNum:
					#print("found entry")
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
	
# Returns the sequence associated with the given miRNA name	
def mirnaSeq(mirnaName):
	for line in mirnaInfo:
		mirnaCmp = line[0]
		mirnaCmp = mirnaCmp.split(" ")
		mirnaCmp = mirnaCmp[0]
		if mirnaCmp == mirnaName:
			return str(line[1])
			
def runMiranda(reprocessLine):
	# Create temp input text files for miRNA and SNP fasta seqs
	# Run miranda on each pair using temp input files
	# Output to temp output text file
	# Delete temp input text files
	# Open output text file, store line, delete output text file
	
	toReturn = None
		
	outputFile = "temp_mirna_input.fasta"
	with open(outputFile, "a") as text_file:
		header = reprocessLine[0]
		sequence = reprocessLine[1]
		
		# Print to file 
		print("{}".format(header), file=text_file)
		print("{}".format(sequence), file=text_file)
		
	outputFile = "temp_snp_input.fasta"
	with open(outputFile, "a") as text_file:
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
		"temp_mirna_input.fasta", 
		"temp_snp_input.fasta", 
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
	
	if toReturn==None:
		print(mirandaText)
	
	# Delete temp input files
	toRun = ["rm", "temp_mirna_input.fasta", "temp_snp_input.fasta"]
	subprocess.run(toRun, check=True)
	
	'''
	print(mirandaTextArray)
	return
	
	for line in f:
		if line[0:2]=='>h':
			toReturn = line.rstrip()
	'''
	
	# Return the score line
	return toReturn