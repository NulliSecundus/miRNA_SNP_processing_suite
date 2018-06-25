import click

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
        refInfo = []
        with open(chromReport) as f:
            for line in f:
                lineSplit = line.split("\t")
                refInfo.append(lineSplit[:13])
        refInfo = refInfo[7:]
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
            for line in f:
                if line[0]==">":
                    sequence = ""
                    header = line
                    header = header.replace('\n', '')
                    headerSplt = header.split("|")
                    rsNum = headerSplt[2].split(" ")
                    rsNum = rsNum[0][2:]

                    unique = (int(refInfo[index][1])==2)
                    gene = (len(refInfo[index][12])>1)

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

                    index += 1
                elif line[0]=="\n":
                    if unique and gene and stdSNP:
                        with open(outputFile, "a") as text_file:
                            sequence = sequence.replace('\n', '')
                            sequence = sequence.replace(' ', '')
                            sequence += '\n'
                            
                            for allele in alleles:
                                if pos>=26:
                                    seqChop = sequence[pos-26:pos-1]
                                #seqChop += " "
                                if allele=='-':
                                    seqChop += ''
                                else:
                                    seqChop += allele
                                    #seqChop += " "
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
