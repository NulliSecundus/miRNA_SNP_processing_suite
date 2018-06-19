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
		with open(procSnpFasta) as f:
			for line in f:
				if line[0]==">":
					header = line
					header = header.replace('\n', '')
					print("header")
					
					snpName = header.split(" ")
					snpName = snpName[0][1:]
					print("snp name")
                    
					headerSplt = header.split("|")
					alleles = headerSplt[8]
					alleles = alleles[8:]
					alleles = alleles.replace('\"','')
					alleles = alleles.split('/')
					alleleNum = len(alleles)
					print("alleles")
					
					temp = [snpName, alleleNum]
					snpInfo.append(temp)
					print("append")
					
					if count == 10:
						print(snpInfo)
						return
					
					count = count+1
				else:
					pass
	except:
		print('Could not parse processed snp fasta file')

def iterateMiranda():
	try: 
		print("success")
	except:
		print("error")