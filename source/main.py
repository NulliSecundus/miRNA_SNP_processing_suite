'''
Created on Aug 22, 2017

@author: JayYang
'''

if __name__ == '__main__':
    pass

container = [0]
sequence = ""
header = ""

# file input name - change as needed
with open("rs_chY.fas") as f:
    for line in f:
        #do something
        # print(line)
        if line[0]==">":
            sequence = ""
            container = [0]
            container[0] = line
        elif line[0]=="\n":
            for i in container:
                if i[0]!=">":
                    sequence += i
            with open("snp_list.fasta", "a") as text_file:
                header = container[0]
                header = header.replace('\n', '')
                headerSplt = header.split("|")
                
                pos = headerSplt[3]
                pos = int(pos[4:])
                
                alleles = headerSplt[8]
                alleles = alleles[8:]
                alleles = alleles.replace('\"','')
                alleles = alleles.split('/')
                
                sequence = sequence.replace('\n', '')
                sequence = sequence.replace(' ', '')
                sequence += '\n'
                #print(sequence[pos-1])
                
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
            comment = line
        else:
            container.append(line)