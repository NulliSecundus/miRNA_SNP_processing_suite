import click

refInfo = []
refGenes = []
index = 0
chromosome = -1
s = False

@click.command()
@click.argument('refflat')
@click.argument('reportfile')
@click.argument('snpfile')
@click.argument('output')
@click.option('--strict', is_flag=True, help='''Strict filtering''')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(refflat, reportfile, snpfile, output, strict, verbose):
	"""
	\b
	Arguments: 
		REFFLAT (refFlat.txt) - INPUT gene annotation file. 
		REPORTFILE (chr_N.txt) - INPUT chromosome report file containing SNP 
			reference information. 
			See ftp://ftp.ncbi.nlm.nih.gov/snp/00readme.txt 
		SNPFILE (rs_chN.fas) - INPUT SNP sequence files by chromosome in fasta 
			format. See ftp://ftp.ncbi.nlm.nih.gov/snp/00readme.txt 
		OUTPUT (proc_chrN.fasta) - processed fasta formatted file 
			containing entries with sequences for each validated SNP. 
	
	Extracts validated SNPs with a unique genomic position located within 
	the 5' or 3' UTR or CDS of any gene. Trims rs_chrN.fas entries to 25 bp
	each side of SNP and creates separate FASTA entries for each allele."""
	global s
	if strict:
		s = True
	try:
		loadRef(refflat)
		loadReport(reportfile)
		procSNP(snpfile, output, verbose)
		print("Success")
	except Exception as error:
		toPrint = "Error: " + str(error)
		print(toPrint)
		return
		
# Function for loading refFlat.txt
def loadRef(refFlat):
	global refGenes
	
	for x in range(0,24):
		refGenes.append([])
	
	print("Loading refFlat.txt")
	
	with open(refFlat) as f:
		for line in f:
			lineSplit = line.split("\t") # Split line by tabs
			geneName = lineSplit[0] # First item is the gene's name
			chr = lineSplit[2] # Chromosome number where gene is located
			
			try:
				if len(chr) > 5:
					chrSplit = chr.split("_")
					chr = chrSplit[0]
				
				chrNum = -1
				chr = chr[3:]
				if chr == "X":
					chrNum = 22
				elif chr == "Y":
					chrNum = 23
				else:
					chrNum = int(chr)-1
				refGenes[chrNum].append(geneName)
			except:
				continue

# Function for loading the chromosome report file 
def loadReport(chromReport):
	global refInfo
	global chromosome
	
	print("Loading SNP reference info")
	
	# Open the chromosome report file to extract reference info 
	lineSplit = None
	with open(chromReport) as f:
		# Skip the first 7 lines
		for x in range(7):
			next(f)
			
		num = 0
		# Starting at line 8
		for line in f:
			lineSplit = line.split("\t") # Split line by tabs
			#lineSplit = lineSplit[:13] # Keep parts 0-12 only 
			tempLine = [
				int(lineSplit[0]), # rs number
				int(lineSplit[1]), # Map (2 == mapped to single location in genome)
				int(lineSplit[2]), # SNP Type (0 == not withdrawn)
				int(lineSplit[3]), # Total number of chromosomes hit
				int(lineSplit[4]), # Total number of contigs hit
				int(lineSplit[5]), # Total number of hits to genome
				lineSplit[12], # Genes at this same position on the chromosome 
				int(lineSplit[16]), # Validated status (0 == no validation information) 
				int(lineSplit[17]), # Genotypes available in dbSNP for this RefSNP
				int(lineSplit[18]), # Linkout available to submitter website for further data on the RefSNP
				lineSplit[25].replace('\n', '') # Global Minor Allele Frequency (GMAF)
			]
			refInfo.append(tempLine)
			'''
			print(tempLine)
			num += 1
			if (num%16==0):
				return
			'''
			
			
	chr = lineSplit[6]
	if chr == "X":
		chromosome = 22
	elif chr == "Y":
		chromosome = 23
	else:
		chromosome = int(chr)-1
	
			
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
				if validateGene(rsNum) and stdSNP and (pos>25):
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
	
def validateGene(rsNum):
	global refGenes
	global chromosome
	
	tempLine = rsSearch(rsNum) # Search for rs num, return info
	if tempLine == None :
		# Terminate program if nothing is found
		print("Error")
		return
	else:
		# Otherwise assign values and continue 
		unique = (tempLine[1]==2)
		gene = str(tempLine[6])
	
	if not unique:
		# SNP must be unique (mapped to single location in chromosome)
		# If not, return False
		return False
	
	if len(gene)<1 :
		# SNP must be mapped to a gene
		# If not, return False
		return False
		
	if s:
		# Strict filtering section
		if (tempLine[2]==0):
			return False
		if (tempLine[3]!=1):
			return False
		if (tempLine[4]!=1):
			return False
		if (tempLine[5]!=1):
			return False
		if (tempLine[7]==0):
			return False
		if (tempLine[8]!=1):
			return False
		if (tempLine[9]!=1):
			return False
		if (len(tempLine[10])<5):
			return False
	
	if gene in refGenes[chromosome]:
		# Mapped gene must be in refGenes list
		return True
		
	# If mapped gene cannot be found in refGenes, return False
	return False