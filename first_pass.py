import click
import subprocess
import secrets
from multiprocessing import Pool

snpList = []
rnaList = []
mirandaList = []
sigID = None
dir = None
sc = None 

@click.command()
@click.argument('snpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.argument('score')
def cli(snpfile, mirnafile, output, score):
	global sc 
	sc = score
	
	try:
		genSig()
		loadsnp(snpfile)
		loadrna(mirnafile)
		iterateMiranda(output)
	except:
		print("Error")
		return
	
# Generates a unique signature ID to be used in file naming
def genSig():
	global sigID
	global dir 
	
	sigID = str(secrets.randbelow(999999999999))
	dir = "temp_" + sigID
	
	toRun = ["mkdir", dir]
	subprocess.run(toRun, check=True)
	
	dir += "/"

# Loads the SNP file into memory and splits it every 3 million entries
# Outputs each split section as a temp file for miranda 	
def loadsnp(snpFile):
	global snpList
	
	print("Loading SNP fasta file")
	
	snpBlock = None
	fileNum = 1
	count = 0
	
	with open(snpFile) as f:
		for line in f:
			if line[0]==">":
				snpBlock = [line.replace('\n', ''), "seq"]
			elif line[0]=="\n": 	
				if snpBlock != None:
					snpList.append(snpBlock)
					count += 1
					if count%3000000==0:
						outputSnp(fileNum)
						fileNum += 1
			else:
				snpBlock[1] = line.replace('\n', '')
	
	# Output any remaining SNP entries
	if count%3000000!=0:
		outputSnp(fileNum)
	
# Loads the miRNA file into memory and splits it every 200 entries 
# Outputs each split section as a temp file for miranda 
def loadrna(mirnaFile):
	global rnaList
	
	print("Loading miRNA fasta file")
	
	rnaBlock = None
	fileNum = 1
	count = 0
	
	with open(mirnaFile) as f:
		for line in f:
			if line[0]==">":
				rnaBlock = [line.replace('\n', ''), "seq"]
			else:
				rnaBlock[1] = line.replace('\n', '')
				rnaList.append(rnaBlock)
				count += 1
				if count%200==0:
					outputRna(fileNum)
					fileNum += 1
				
	# Output any remaining miRNA entries
	if count%200!=0:
		outputRna(fileNum)
	
def iterateMiranda(out):
	print("Running miranda on SNP part 1")
	iterate = []
	for i in range(15):
		temp = [str(i), " test", " result"]
		iterate.append(temp)
	with Pool() as p:
		p.map(test, iterate)
		
def test(x):
	result = x[0] + x[1] + x[2]
	print(result)
	
# Utility function for printing current snpList to file
def outputSnp(n):
	global snpList
	
	tempSnpFileName = dir + sigID + "_snp_" + str(n) + ".fasta"
	with open(tempSnpFileName, "w") as text_file:
		for entry in snpList:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	snpList = []
	
# Utility function for printing current rnaList to file
def outputRna(n):
	global rnaList
	
	tempRnaFileName = dir + sigID + "_mirna_" + str(n) + ".fasta"
	with open(tempRnaFileName, "w") as text_file:
		for entry in rnaList:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	rnaList = []