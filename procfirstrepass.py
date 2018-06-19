import click

@click.command()
@click.argument('mirandafile')
@click.argument('procsnpfile')
@click.argument('mirna')
@click.argument('output')
def cli(mirandafile, procsnpfile, mirna, output):
	try:
		# Run processInput on original miranda output
		loadsnp(procsnpfile)
	except:
		print("Failed to parse snp file")
	try:
		# Run iterateMiranda on reprocess list
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

def iterateMiranda():
	try: 
		print("success")
	except:
		print("error")