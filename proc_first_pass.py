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
	snpID = "default" # SNP ID
	
	with open(mirandaOut) as f, open(outputFile, "a") as o:
		for line in f:
			if line[0]==">":
				# Data line 
				lineSplit = line.split("\t")
				snpInfo = lineSplit[1].split("|")
				cmpID = snpInfo[2]
				if snpID != cmpID:
					if snp != None:
						# Count number of SNP alleles in output 
						# Compare to total number of SNP alleles 
						snpNum = len(snp)
						
						snpTot = int(snpInfo[3])
						
						if snpNum < snpTot:
							# Print to file 
							for entry in snp:
								print("{}".format(entry), file=o)
						
					# Start new SNP container 
					snp = [line.replace("\n", "")]
					snpID = cmpID
					
				else: 
					snp.append(line.replace("\n", ""))
			elif line[0]=="#":
				# Comment line, append to output file
				print("{}".format(line.replace("\n", "")), file=o)