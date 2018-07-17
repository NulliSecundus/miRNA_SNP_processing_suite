import click
refInfo = []

@click.command()
@click.argument('reportfile')
@click.argument('snpfile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console''')
def cli(reportfile, snpfile, output, verbose):
	try:
		processInput(reportfile, snpfile, output, verbose)
	except:
		print('Error processing the chromosome fasta file.')

if __name__ == '__main__':
	pass

def processInput(chromReport, snpFasta, outputFile, v):
	try:
		with open(chromReport) as f:
			for x in range(7):
				next(f)
			for line in f:
				lineSplit = line.split("\t")
				lineSplit = lineSplit[:13]
				tempLine = [
					int(lineSplit[0]),
					int(lineSplit[1]),
					lineSplit[12]
				]
				refInfo.append(tempLine)
		sequence = ""
		header = ""
		index = 0
		unique = False
		gene = False
		stdSNP = False
	except:
		print('Could not read chromosome report file')
    
	try:
		count = 0
		with open(snpFasta) as f:
			with open(outputFile, "a") as text_file:
				for line in f:
					if line[0]==">":
						sequence = ""
						header = line
						header = header.replace('\n', '')
						headerSplt = header.split("|")
						rs = headerSplt[2].split(" ")
						rsNum = int(rs[0][2:])
						print(rsNum)
						print(refInfo[0])

						tempLine = rsSearch(rsNum)
						print("found rs")
						if tempLine == None :
							unique = False
							gene = False
						else:
							unique = (tempLine[1]==2)
							gene = (len(tempLine[2])>1)

						print("assigned truth values")
						pos = headerSplt[3]
						pos = int(pos[4:])

						snpClass = headerSplt[7]
						snpClass = snpClass[6:]
						stdSNP = (snpClass == '1')
							
						alleles = headerSplt[8]
						alleles = alleles[8:]
						alleles = alleles.replace('\"','')
						alleles = alleles.split('/')
						alleleNum = len(alleles)
						'''
						if int(rsNum) == 62247966:
							print("found")
							print(unique)
							print(gene)
							print(stdSNP)
							print(pos>25)
							print(refInfo[index])
							return
						'''
						
						index += 1
					elif line[0]=="\n":
						if unique and gene and stdSNP and (pos>25):
						#if stdSNP and (pos>25):
							sequence = sequence.replace('\n', '')
							sequence = sequence.replace(' ', '')
							sequence += '\n'
							
							if len(sequence) < 50:
								continue
								
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
						#Do nothing, it's a comment line
						pass
					else:
						sequence += line
	except:
		print('Could not read snp fasta file')
		
def rsSearch(rsNumber):
	print("rs search")
	
	for item in refInfo:
		print("comparing: ")
		print(item[0])
		print(rsNumber)
		if item[0]==rsNumber:
			return item
			
	return None