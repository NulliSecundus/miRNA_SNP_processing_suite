import click

refInfo = []
index = 0

@click.command()
@click.argument('reportfile')
@click.argument('snpfile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(reportfile, snpfile, output, verbose):
	try:
		processInput(reportfile, snpfile, output, verbose)
		print("Success")
	except:
		print("Error")

def processInput(chromReport, snpFasta, outputFile, v):
	global refInfo
	
	print("Processing SNPs in chromosome")
	
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
					gene = (len(tempLine[2])>1)

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
				if unique and gene and stdSNP and (pos>25):
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