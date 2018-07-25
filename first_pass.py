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
				if snpBlock != None:
					snpList.append(snpBlock)
					count += 1
					if count%3000000==0:
						outputSnp(fileNum)
						fileNum += 1
				snpBlock = [line.replace('\n', ''), "seq"]
			elif line[0]=="\n": 	
				pass
			else:
				snpBlock[1] = line.replace('\n', '')
				
	outputSnp(fileNum)
	
# Loads the miRNA file into memory and splits it every 200 entries 
# Outputs each split section as a temp file for miranda 
def loadrna(mirnaFile):
	print("Loading miRNA fasta file")
	
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
	tempSnpFileName = sigID + "_snp_" + str(n) + ".fasta"
	print(tempSnpFileName)
	snpList = []