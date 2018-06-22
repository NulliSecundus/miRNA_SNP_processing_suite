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
		print("Loading SNP fasta file (might take a few minutes)")
		
		count = 0
		header = ""
		snpName = ""
		allele = ""
		alleleNum = 0
		alleleIndex = 0
		sequence = ""
		
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
					temp = [snpName, alleleNum, allele, sequence]
					snpInfo.append(temp)
					'''
					if count%1000000 == 0:
						print(count)
						
						for line in snpInfo:
							print(line)
						
						return
						'''
					
					count = count+1
					
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
		
		tempCont = []
		count = 0
		prevLine = ""
		
		with open(mirandaFile) as f:
			for line in f:
				if line[0]==">":
					lineEdit = line.split("\t")
					mirnaName = lineEdit[0]
					refName = lineEdit[1]
					
					if refName == prevLine:
						tempCont.append(refName)
						prevLine = refName
					else:
						alleleCount = len(tempCont)
						#print(tempCont)
						
						if alleleCount > 0:
							checkAlleleCount(tempCont[0], alleleCount, mirnaName)
						
						tempCont = [refName]
						prevLine = refName
					
					count = count+1
				else:
					pass
					
				if count%10000 == 30:
					print(count)
					
					for entry in reprocessList:
						print(entry)
					
					return
					
	except:
		print('Could not build reprocess list from miranda file')

def iterateMiranda():
	try: 
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
		
def checkAlleleCount(name, num, mirna):
	#print(name + " " + str(num))
	startSig = False 
	temp = [mirna, mirnaSeq(mirna)]
	
	for line in snpInfo:
		strcmp = line[0]
		if strcmp == name:
			if (line[1] - num) > 0:
				startSig = True 
				#print("re-process " + name)
				refSeqName = line[0] + "|" + str(line[1]) + "|" + line[2]
				# Add in other alleles 
				temp += [refSeqName, line[3]]
					
		else: 
			if startSig:
				reprocessList.append(temp)
				return
			
def mirnaSeq(mirnaName):
	#print("Start mirnaSeq")
	for line in mirnaInfo:
		mirnaCmp = line[0]
		mirnaCmp = mirnaCmp.split(" ")
		mirnaCmp = mirnaCmp[0]
		#print(mirnaCmp)
		if mirnaCmp == mirnaName:
			#print("Found miRNA")
			return str(line[1])