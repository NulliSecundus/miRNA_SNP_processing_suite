import click

refInfo = []
refGenes = []
index = 0

@click.command()
@click.argument('refflat')
@click.argument('reportfile')
@click.argument('snpfile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(refflat, reportfile, snpfile, output, verbose):
	"""
	\b
	Arguments: 
		REFFLAT (refFlat.txt) - INPUT gene annotation file. 
		REPORTFILE (chr_N.txt) - INPUT chromosome report file containing SNP 
			reference information. 
			See ftp://ftp.ncbi.nlm.nih.gov/snp/00readme.txt 
		SNPFILE (rs_chrN.fas) - INPUT SNP sequence files by chromosome in fasta 
			format. See ftp://ftp.ncbi.nlm.nih.gov/snp/00readme.txt 
		OUTPUT (proc_chrN.fasta) - processed fasta formatted file 
			containing entries with sequences for each validated SNP. 
	
	Extracts validated SNPs with a unique genomic position located within 
	the 5' or 3' UTR or CDS of any gene. Trims rs_chrN.fas entries to 25 bp
	each side of SNP and creates separate FASTA entries for each allele."""
	try:
		loadRef(refflat)
		print(refGenes[0:5])
		print(len(refGenes))
		loadReport(reportfile)
		print(refInfo[0:5])
		procSNP(snpfile, output, verbose)
		print("Success")
	except:
		print("Error")
		
# Function for loading refFlat.txt
def loadRef(refFlat):
	global refGenes
	
	print("Loading refFlat.txt")
	
	with open(refFlat) as f:
		for line in f:
			lineSplit = line.split("\t") # Split line by tabs
			geneName = lineSplit[0] # First item is the gene's name
			refGenes.append(geneName)

# Function for loading the chromosome report file 
def loadReport(chromReport):
	global refInfo
	
	print("Loading SNP reference info")
	
	# Open the chromosome report file to extract reference info 
	with open(chromReport) as f:
		# Skip the first 7 lines
		for x in range(7):
			next(f)
			
		# Starting at line 8
		for line in f:
			lineSplit = line.split("\t") # Split line by tabs
			lineSplit = lineSplit[:13] # Keep parts 0-12 only 
			tempLine = [
				int(lineSplit[0]), # rs number
				int(lineSplit[1]), # Map (2 == mapped to single location in genome)
				lineSplit[12] # Genes at this same position on the chromosome 
			]
			refInfo.append(tempLine)
			
# Function for processing the SNP fasta file	
def procSNP(snpFasta, outputFile, v):
	print("Building processed SNP list")
	
	# Initialize vars 
	sequence = ""
	header = ""
	unique = False
	gene = False
	stdSNP = False
    
	count = 0
	# Open chromosome SNP fasta file and output file for printing 
	with open(snpFasta) as f, open(outputFile, "w") as text_file:
		for line in f:
			if line[0]==">":
				# Header line 
				# Format and split by "|"
				sequence = ""
				header = line
				header = header.replace('\n', '')
				headerSplt = header.split("|")
				
				rs = headerSplt[2].split(" ")
				rsNum = int(rs[0][2:]) # Extract rs number (INT)

				tempLine = rsSearch(rsNum) # Search for rs num, return info
				if tempLine == None :
					# Terminate program if nothing is found
					print("Error")
					return
				else:
					# Otherwise assign values and continue 
					unique = (tempLine[1]==2)
					gene = str(tempLine[2])
					g = (len(gene)>1)

				pos = headerSplt[3]
				pos = int(pos[4:]) # Position of SNP in sequence 

				snpClass = headerSplt[7]
				snpClass = snpClass[6:]
				stdSNP = (snpClass == '1') # SNP class 1 indicates a standard SNP 
					
				alleles = headerSplt[8]
				alleles = alleles[8:]
				alleles = alleles.replace('\"','')
				alleles = alleles.split('/') # List of SNP alleles 
				alleleNum = len(alleles) # Total number of SNP alleles 
				
			elif line[0]=="\n":
				# Newline indicates moving to next SNP 
				if validateGene(unique, gene) and stdSNP and (pos>25):
					# If prev SNP satisfies criteria
					# Begin formatting for printing 
					
					sequence = sequence.replace('\n', '') # Remove embedded newlines
					sequence = sequence.replace(' ', '') # Remove spaces
					sequence += '\n' # Add a single newline to end of sequence 
					
					# Skip if lenth of sequence is less than 50
					if len(sequence) < 50:
						continue
						
					# Create a separate fasta entry for each SNP allele 
					for allele in alleles:
						seqChop = sequence[pos-26:pos-1]
						if allele=='-':
							seqChop += ''
						else:
							seqChop += allele
						seqChop += sequence[pos:pos+25]
						seqChop += "\n" 
						
						# Modify header with allele info
						label = header.split(" ")
						label = label[0]
						label += "|" + str(alleleNum)
						label += "|" + allele
						
						# Print each header/sequence to file 
						print("{}".format(label), file=text_file)
						print("{}".format(seqChop), file=text_file)
						
						# Increment counter
						count += 1
						if v:
							if count%10000==0:
								print(count)
			elif line[0]=="#":
				# Do nothing, it's a comment line
				pass
			else:
				# Add line to the sequence entry 
				sequence += line
	
# Utility function for searching refInfo for the specified rs number 
def rsSearch(rsNumber):
	global index
	global refInfo
	
	for n in range(index, len(refInfo)):
		if refInfo[n][0]==rsNumber:
			index = n 
			return refInfo[n] # Return whole dataline 
	return None
	
def validateGene(unique, gene):
	global refGenes
	
	if not unique:
		# SNP must be unique (mapped to single location in chromosome)
		# If not, return False
		return False
	
	if len(gene)<1 :
		# SNP must be mapped to a gene
		# If not, return False
		return False
	
	if gene in refGenes:
		# Mapped gene must be in refGenes list
		return True
		
	# If mapped gene cannot be found in refGenes, return False
	return False