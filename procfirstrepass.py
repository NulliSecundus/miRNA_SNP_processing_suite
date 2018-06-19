import click

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('output')
def cli(mirandafile, procsnpfile, output):
	try:
		# Run processInput on original miranda output
		processInput(mirandafile, procsnpfile)
	except:
		print("Failed to parse miranda file")
	try:
		# Run iterateMiranda on reprocess list
		iterateMiranda()
	except:
		printo("Failed to re-process with miranda")

if __name__ == '__main__':
	pass

reprocessList = []
	
def processInput(mirandaFile, procSnpFasta):
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
					
					if count%10000 == 0:
						print(count)
						print(snpInfo[0])
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

def iterateMiranda():
	try: 
		print("success")
	except:
		print("error")