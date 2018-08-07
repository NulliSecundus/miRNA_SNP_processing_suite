import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
def cli(mirandafile, output):
	try:
		parsemir(mirandafile, output)
		print("Success")
	except:
		print("Error")

def parsemir(mirandaOut, outputFile):
	print("Processing first pass miranda output file")
	
	snp = None # SNP container 
	snpInfo = None # SNP ref info
	snpID = "default" # SNP ID
	count = 1
	
	with open(mirandaOut) as f, open(outputFile, "w") as o:
		for line in f:
			if line[0]==">":
				# Data line 
				print(line)
				
				lineSplit = line.split("\t") # Split data line by tabs
				
				print(lineSplit)
				
				snpInfo = lineSplit[1].split("|") # Split the SNP entry ID
				cmpID = snpInfo[2] # Use the rs number from the SNP entry ID 
				
				toPrint = "Compare ID: " + cmpID
				print(toPrint)
				toPrint = "SNP ID: " + snpID
				print(toPrint)
				print("\n")
				
				if snpID != cmpID:
					if snp != None:
						# For every case except the start 
						snpNum = len(snp) # Count number of SNP alleles in output 
						snpTot = int(snpInfo[3]) # Get total num of alleles 
						
						# Compare to total number of SNP alleles 
						if snpNum < snpTot: 
							# Print to file 
							for entry in snp:
								print("{}".format(entry), file=o)
						
					# Start new SNP container 
					snp = [line.replace("\n", "")]
					snpID = cmpID # Set snpID to the new rs number
					
				else: 
					snp.append(line.replace("\n", ""))
				
				if count%20==0:
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