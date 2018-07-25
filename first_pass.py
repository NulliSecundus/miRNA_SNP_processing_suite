import click
import subprocess
import secrets
from multiprocessing import Pool

snpList = []
rnaList = []
mirandaList = []
sigID = None

@click.command()
@click.argument('snpfile')
@click.argument('mirnafile')
@click.argument('output')
@click.argument('score')
def cli(snpfile, mirnafile, output, score):
	try:
		genSig()
		loadsnp(snpfile)
		loadrna(mirnafile)
		iterateMiranda(output, score)
	except:
		print("Error")
		return
	
def genSig():
	global sigID
	sigID = str(secrets.randbelow(999999999999))
	
def loadsnp(snpFile):
	print("Loading SNP fasta file")
	
def loadrna(mirnaFile):
	print("Loading miRNA fasta file")
	
def iterateMiranda(out, sc):
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
	for n in range(999999999):
		s = str(n)