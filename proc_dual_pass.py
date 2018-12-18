import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
@click.option('--debug', is_flag=True, help='Show debug messages')
def cli(mirandafile, output, debug):

	"""
	\b
	Arguments:
		MIRANDAFILE - INPUT processed final output from dualpass
		OUTPUT - combined output sqlite file 
	
	Combines the output from a miranda dual pass into a single sqlite file."""

	try:
		parsemir(mirandafile, output)
		print("Success")
	except:
		print("Error")

def parsemir(mirandaOut, outputFile):
	snp = None 
	snpInfo = None 
	snpID = "default" 
	snpTot = 0 
	count = 1
	
	with open(mirandaOut) as f, open(outputFile, "w") as o:
		for line in f:
			if line[0]==">":
				lineSplit = line.split("\t")
				snpInfo = lineSplit[1].split("|") 
				cmpID = snpInfo[2] 
				if snpID != cmpID:
					if snp != None:
						snpNum = len(snp)
						if snpNum < snpTot: 
							for entry in snp:
								print("{}".format(entry), file=o)
					snp = [line.replace("\n", "")]
					snpID = cmpID
					snpTot = int(snpInfo[3])
				else: 
					snp.append(line.replace("\n", ""))
			elif line[0]=="#":
				if line[0:5]=="# End":
					snpNum = len(snp)
					snpTot = int(snpInfo[3])
					
					if snpNum < snpTot:
						for entry in snp:
							print("{}".format(entry), file=o)
				print("{}".format(line.replace("\n", "")), file=o)