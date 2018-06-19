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
	try:
		# Run loadrna
		loadrna(mirna)
	except:
		print("Failed to load miRNA file")
	try:
		# Run buildReproc on original miranda output
		buildReproc(mirandafile)
	except:
		print("Failed to parse miranda file")
	try:
		# Run iterateMiranda to reprocess list
		iterateMiranda()
	except:
		printo("Failed to re-process with miranda")

if __name__ == '__main__':
	pass

reprocessList = []
	
def loadsnp(procSnpFasta):
	snpInfo = []
	
	try:
		count = 0
		header = ""
		snpName = ""
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
					
					snpName = snpName + "|" + alleles[alleleIndex]
					
					alleleIndex = alleleIndex + 1
					
				elif line[0]=="\n":
					# End of sequence
					temp = [snpName, alleleNum, sequence]
					snpInfo.append(temp)
					
					if count%10000 == 10:
						print(count)
						print(snpInfo[0])
						print(snpInfo[1])
						print(snpInfo[2])
						print(snpInfo[3])
						print(snpInfo[4])
						print(snpInfo[5])
						print(snpInfo[6])
						return
					
					count = count+1
					
				elif line[0]=="#":
					#Do nothing, it's a comment line
					pass
				else:
					# Sequence line 
					sequence = line

	except:
		print('Could not parse processed snp fasta file')
	
def loadrna(miRNA):
	mirnaInfo = []
	
	try:
		count = 0
		header = ""
		sequence = ""
		
		with open(miRNA) as f:
			for line in f:
				if line[0]==">":
					header = line
					
				else:
					# Sequence line 
					sequence = line
					
					temp = [header, sequence]
					mirnaInfo.append(temp)
					
					if count%10000 == 10:
						print(count)
						print(mirnaInfo[0])
						print(mirnaInfo[1])
						print(mirnaInfo[2])
						print(mirnaInfo[3])
						print(mirnaInfo[4])
						print(mirnaInfo[5])
						print(mirnaInfo[6])
						return
					
					count = count+1

	except:
		print('Could not parse miRNA file')

def buildReproc(mirandaFile):
	try: 
		print("success")
		
		# For each line in cleaned miranda output
		# Populate temp container with same-name entries
		# When new name appears, count num in temp container
		# Search snpInfo for the entry, compare to temp container num
		# If (available - output) > 0, then add SNP-miRNA pair to reprocess list
		
	except:
		print("error")

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