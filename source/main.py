'''
Created on Aug 22, 2017

@author: Brendan Drackley
'''

if __name__ == '__main__':
    pass

chromReport = "chr_1.txt"
snpFasta = "rs_ch1.fas"

refInfo = []
with open(chromReport) as f:
    for line in f:
        lineSplit = line.split("\t")
        refInfo.append(lineSplit[:2])
refInfo = refInfo[7:]
sequence = ""
header = ""
gate = False
index = 0

# file input name - change as needed
with open(snpFasta) as f:
    for line in f:
        if line[0]==">":
            sequence = ""
            header = line
            header = header.replace('\n', '')
            headerSplt = header.split("|")
            rsNum = headerSplt[2].split(" ")
            rsNum = rsNum[0][2:]
            gate = (int(refInfo[index][1])==2)
                
            pos = headerSplt[3]
            pos = int(pos[4:])
                
            alleles = headerSplt[8]
            alleles = alleles[8:]
            alleles = alleles.replace('\"','')
            alleles = alleles.split('/')
            index += 1
        elif line[0]=="\n":
            if gate:
                with open("snp_list_chrY.fasta", "a") as text_file:
                    sequence = sequence.replace('\n', '')
                    sequence = sequence.replace(' ', '')
                    sequence += '\n'
                    
                    for allele in alleles:
                        if pos>=26:
                            seqChop = sequence[pos-26:pos-1]
                        seqChop += " "
                        if allele=='-':
                            seqChop += ''
                        else:
                            seqChop += allele
                            seqChop += " "
                        seqChop += sequence[pos:pos+25]
                        seqChop += "\n"
                        print("{}".format(header), file=text_file)
                        print("{}".format(seqChop), file=text_file)
        elif line[0]=="#":
            #Do nothing, it's a comment
            pass
        else:
            sequence += line