import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
@click.option('-stop', default=20, help='''
Limits processing to the specified number of SNP entries\n
Parameter only used if debug is flagged\n
Default: 20
''')
@click.option('--debug', is_flag=True, help='Show debug messages')
def cli(mirandafile, output, stop, debug):
	try:
		parsemir(mirandafile, output, stop, debug)
		print("Success")
	except:
		print("Error")

def parsemir(mirandaOut, outputFile, stop, d):
	print("Processing first pass miranda output file")
	
	snp = None # SNP container 
	snpInfo = None # SNP ref info
	snpID = "default" # SNP ID
	snpTot = 0 # Total number of SNP alleles
	count = 1
	
	with open(mirandaOut) as f, open(outputFile, "w") as o:
		for line in f:
			if line[0]==">":
				# Data line 
				if d:
					toPrint = "Count number: " + str(count)
					print(toPrint)
					print(line)
				
				lineSplit = line.split("\t") # Split data line by tabs
				
				if d:
					print(lineSplit)
				
				snpInfo = lineSplit[1].split("|") # Split the SNP entry ID
				cmpID = snpInfo[2] # Use the rs number from the SNP entry ID 
				
				if d:
					toPrint = "Compare ID: " + cmpID
					print(toPrint)
					toPrint = "SNP ID: " + snpID
					print(toPrint)
				
				if snpID != cmpID:
					if snp != None:
						# For every case except the start 
						snpNum = len(snp) # Count number of SNP alleles in output 
						
						if d:
							toPrint = "SNP Alleles: " + str(snpNum) + "\nTotal Alleles: " + str(snpTot)
							print(toPrint)
						
						# Compare to total number of SNP alleles 
						if snpNum < snpTot: 
							# Print to file 
							if d:
								print("Printing previous SNPs:")
							
							for entry in snp:
								print("{}".format(entry), file=o)
								
								if d:
									print(entry)
						
					# Start new SNP container and info 
					snp = [line.replace("\n", "")]
					snpID = cmpID # Set snpID to the new rs number
					snpTot = int(snpInfo[3]) # Get new total num of alleles 
					
				else: 
					if d:
						print("Appending SNP line to container")
					
					snp.append(line.replace("\n", ""))
				
				if d:
					print("\nEnd of line")
					print("\n")
					print("\n")
					
					if count%stop==0:
						return
					
					count += 1
					
			elif line[0]=="#":
				if line[0:5]=="# End":
					# Compare SNP alleles for the last SNP container 
					snpNum = len(snp)
					snpTot = int(snpInfo[3])
					
					if snpNum < snpTot:
						# Print to file 
						for entry in snp:
							print("{}".format(entry), file=o)
							
				# Comment line, append to output file
				print("{}".format(line.replace("\n", "")), file=o)