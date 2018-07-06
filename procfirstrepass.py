import click

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('mirna')
@click.argument('output')
def cli(mirandafile, procsnpfile, mirna, output):
	try:
		# Run loadsnp
		loadsnp(procsnpfile)
		print(len(snpInfo))
		print(len(snpInfo[0]))
		print(snpInfo[0][1])
		print(snpInfo[0][0])
		print(snpInfo[1][0])
		print(snpInfo[2][0])
		
	except:
		print("Failed to load snp file")
		return
	try:
		# Run loadrna
		loadrna(mirna)
	except:
		print("Failed to load miRNA file")
		return
	try:
		# Run buildReprocList on original miranda output
		buildReprocList(mirandafile)
	except:
		print("Failed to parse miranda file")
		return
	try:
		# Run addSequences to add sequences to reprocess list 
		addSequences()
	except:
		print("Failed to add sequences to re-process list")
		return
	try:
		# Run iterateMiranda to reprocess list
		iterateMiranda()
	except:
		printo("Failed to re-process with miranda")
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
					temp = [snpName, alleleNum, allele, sequence, int(rsNum)]
					snpSubInfo.append(temp)
					'''
					if count%1000000 == 0:
						print(count)
						
						for line in snpInfo:
							print(line)
						
						return
						'''
					count = count+1
					
					if count%10000 == 0:
						rsEnd = int(rsNum)
						headerLine = [rsStart, rsEnd]
						snpSubInfo.insert(0, headerLine)
						snpInfo.append(snpSubInfo)
						
						# Reset sub-list and rsStart number
						snpSubInfo = []
						rsStart = rsNum
					
				elif line[0]=="#":
					#Do nothing, it's a comment line
					pass
				else:
					# Sequence line 
					sequence = line
					sequence = sequence.replace('\n', '')

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
					
		# TODO: reorganize reprocessList into twenty subsections
		
	except:
		print('Could not build reprocess list from miranda file')

def addSequences():
	print('Loading sequences into reprocess list')
	print(len(reprocessList))
	mirnaName = ""
	snpName = ""
	count = 0
	
	for line in reprocessList:
		mirnaName = line[0]
		snpName = line[1]
		line.insert(1, mirnaSeq(mirnaName))
		temp = [snpSeq(snpName)]
		line.extend(temp)
		count += 1
		
		if count%50==0:
			print(count)
			print(reprocessList[count-1])
		
def iterateMiranda():
	try: 
		print(reprocessList[0])
		print(len(reprocessList))
		print("success")
		
		# Iterate through list of SNP-miRNA pairs that need to be reprocessed 
		# Create input text files for SNP and miRNA fasta seqs
		# Run miranda on each pair
		# Output to output text file
		# Delete input text files
		# Open output text file, store line, delete output text file
		# Add line to condensed final output file

	except:
		print("error")
		
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

	for line in snpInfo:
		rsStart = line[0][0]
		rsEnd = line[0][1]
		if (rsNum > rsStart) and (rsNum >= rsEnd):
			for entry in line[1:]:
				'''
				snpCmp = entry[0]
				cmpNameSplit = snpCmp.split("|")
				cmpRsText = cmpNameSplit[2]
				cmpRsNum = int(cmpRsText[2:])
				'''
				cmpRsNum = entry[4]
				if cmpRsNum == rsNum:
					return str(entry[3])
					
	print("Failed to locate SNP Sequence")
	
# Returns the sequence associated with the given miRNA name	
def mirnaSeq(mirnaName):
	for line in mirnaInfo:
		mirnaCmp = line[0]
		mirnaCmp = mirnaCmp.split(" ")
		mirnaCmp = mirnaCmp[0]
		if mirnaCmp == mirnaName:
			return str(line[1])