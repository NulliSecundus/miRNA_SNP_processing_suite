import click
import subprocess
import secrets
from multiprocessing import Pool
import time;

snpList = []
snpFileList = []
rnaList = []
rnaFileList = []
mirandaList = [] # Ordered parameters to iterate miranda
outputFileList = []
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
	sc = float(score)
	
	try:
		genSig(snpfile, mirnafile, output)
		loadsnp(snpfile)
		loadrna(mirnafile)
		iterateMiranda(output)
	except:
		print("Error")
		return
	
# Generates a unique signature ID to be used in file naming
def genSig(snpFile, mirnaFile, out):
	global sigID
	global dir 
	
	sigID = str(secrets.randbelow(999999999999))
	dir = "temp_" + sigID
	
	toRun = ["mkdir", dir]
	subprocess.run(toRun, check=True)
	
	dir += "/"
	
	# Create output file and add timestamp
	with open(out, 'w') as o:
		localtime = "# Start: " + time.asctime(time.localtime(time.time()))
		print("{}".format(localtime), file=o)
	
	# Create README info file 
	infoFile = dir + "README" + ".txt"
	with open(infoFile, "w") as text_file:
		header = "Input Parameters: "
		
		# Print to file 
		print("{}".format(header), file=text_file)
		print("{}".format(snpFile), file=text_file)
		print("{}".format(mirnaFile), file=text_file)
		print("{}".format(out), file=text_file)
		print("{}".format(sc), file=text_file)

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
					if count%2000000==0:
						outputSnp(fileNum)
						fileNum += 1
			else:
				snpBlock[1] = line.replace('\n', '')
	
	# Output any remaining SNP entries
	if count%2000000!=0:
		outputSnp(fileNum)
	
# Loads the miRNA file into memory and splits it every 200 entries 
# Outputs each split section as a temp file for miranda 
def loadrna(mirnaFile):
	global rnaList
	global outputFileList
	
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
				if count%80==0:
					outputRna(fileNum)
					fileNum += 1
				
	# Output any remaining miRNA entries
	if count%80!=0:
		outputRna(fileNum)
		
	# Populate list of output files based on miRNA input files
	for n in range(fileNum):
		outName = dir + sigID + "_out_" + str(n+1) + ".txt"
		outputFileList.append(outName)
	
def iterateMiranda(out):
	# Setup lists of files to iterate miranda over 
	for entry in snpFileList:
		tempMiranda = []
		for line in rnaFileList:
			tempMiranda.append([entry, line])
		mirandaList.append(tempMiranda)
	
	# Iteratively run miranda
	num = 1
	for entry in mirandaList:
		toPrint = "Running miranda on SNP part " + str(num)
		print(toPrint)
		with Pool() as p:
			p.map(runMiranda, entry)
		parseMiranda(num)
		appendOutput(out)
		num += 1
		
	# Add ending timestamp to output file 
	with open(out, 'a') as o:
		localtime = "# End: " + time.asctime(time.localtime(time.time()))
		print("{}".format(localtime), file=o)
	
# Utility function for printing current snpList to file
def outputSnp(n):
	global snpList
	global sigID
	global dir
	
	tempSnpFileName = dir + sigID + "_snp_" + str(n) + ".fasta"
	with open(tempSnpFileName, "w") as text_file:
		for entry in snpList:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	snpFileList.append(tempSnpFileName)
	snpList = []
	
# Utility function for printing current rnaList to file
def outputRna(n):
	global rnaList
	global sigID
	global dir
	
	tempRnaFileName = dir + sigID + "_mirna_" + str(n) + ".fasta"
	with open(tempRnaFileName, "w") as text_file:
		for entry in rnaList:
			header = entry[0]
			sequence = entry[1] + "\n"
			
			# Print to file 
			print("{}".format(header), file=text_file)
			print("{}".format(sequence), file=text_file)
	
	rnaFileList.append(tempRnaFileName)
	rnaList = []
	
# Utility function for handling each iteration of miranda 
def runMiranda(x):
	global dir
	global sigID
	
	tempSnpFile = x[0]
	tempRnaFile = x[1]
	
	# Get the file number for output naming
	textArray = tempRnaFile.replace('.fasta', '').split("_")
	fileNum = textArray[3]
	
	# Assign name of temp output file 
	tempOut = dir + sigID + "_out_" + fileNum + ".txt"
	
	# Run miranda
	toRun = [
		"miranda", 
		tempRnaFile, 
		tempSnpFile, 
		"-sc",
		str(sc),
		"-noenergy",
		"-quiet",
		"-keyval",
		"-out",
		tempOut
	]
	subprocess.run(toRun, check=True)

# Utility function for parsing temp miranda outputs into final format 	
def parseMiranda(n):
	toPrint = "Parsing part " + str(n) + " miranda outputs"
	print(toPrint)
	
	for entry in outputFileList:
		output_container = [] # container for lines of final output
		container = [] # temporary container for comparison within a group
		
		with open(entry) as f:
			
			count = 0
			sumcount = 0
			
			for line in f:
				# for each line in the miranda output file
				if line[0:2]=='>h':
					# if line is a data line
					container.append(line.rstrip())
							
				elif line[0:2]=='>>':
					# if line is a summary line
					# end case for a grouping
					# copy top score from from container
					sumcount += 1
					topLine = ""
					
					if len(container)==1 :
						# if the group only has one data line 
						topLine = str(container[0])
					else:
						# if multiple data lines in grouping
						# determine top score line
						# send top score line to output_container
						topScore = -1
						
						for dataline in container:
							# for each data line in the grouping
							splitdline = dataline.split("\t")
							compare = float(splitdline[2])
							if compare > topScore :
								topScore = compare
								topLine = dataline
					
					# append topLine to output_container
					output_container.append(topLine)
						
					# reset the temporary container
					container = [] 
		
		cmpOut = entry.replace('.txt', '_cmp.txt')
		with open(cmpOut, 'w') as o:	
			for line in output_container:
				# write each topLine to the output file 
				print('{}'.format(line), file=o)

		# Delete temp input files
		toRun = ["rm", entry]
		subprocess.run(toRun, check=True)
		
# Utility function for appending each compressed temp output to final
def appendOutput(out):
	for entry in outputFileList:
		cmpOut = entry.replace('.txt', '_cmp.txt')
		with open(cmpOut) as f:
			with open(out, 'a') as o:
				for line in f:
					print('{}'.format(line.replace('\n', '')), file=o)
		
		# Delete temp input files
		toRun = ["rm", cmpOut]
		subprocess.run(toRun, check=True)
		
	# Delete all temp files and folder
	toRun = ["rm", "-rf", dir.replace('/', '')]
	subprocess.run(toRun, check=True)
		