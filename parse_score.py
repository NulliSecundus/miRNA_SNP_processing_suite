import click
import numpy as np

@click.command()
@click.argument('mirandafile')
@click.option('-ut', default=80.0, help='The upper threshold in which the top percentile is calculated\nDefault 80.0')
@click.option('-lt', default=20.0, help='The lower threshold in which the bottom percentile is calculated\nDefault 20.0')
@click.option('-out', default='off', help='Print the parsed scores to the specified output file')
@click.option('--so', is_flag=True, help='Scores only option (for large files). Will only output scores without stats')
@click.option('--scorefile', is_flag=True, help='Flag if a pre-processed score file is provided instead of miranda output')
@click.option('--parsefile', is_flag=True, help='Flag if a parsed miranda output file is provided as input')
def cli(mirandafile, ut, lt, out, so, scorefile, parsefile):
	try:
		if scorefile:
			readScore(mirandafile, ut, lt)
		elif parsefile:
			parseScoreNew(mirandafile, out, so, ut, lt)
		else:
			parseScore(mirandafile, out, ut, lt)
		print('Success')
	except:
		print('Error')

def parseScore(mirandaOut, outputFile, ut, lt):
	scorelist = []
	with open(mirandaOut) as f:
		print('Processing miranda output file')
		for line in f:
			if line[0:2]=='>>':
				lineSplit = line.split('\t')
				scorelist.append(float(lineSplit[4]))
	if outputFile!="off":
		try:
			with open(outputFile, "a") as text_file:
				print('Writing scores to output file')
				for score in scorelist:
					print('{}'.format(score), file=text_file)
		except:
			print('Could not print to output file')
	up = np.percentile(scorelist, ut)  
	print('Top ', ut, ' percentile is: ', up)
	lp = np.percentile(scorelist, lt)  
	print('Bottom ', lt, ' percentile is: ', lp)
	print('Score Parse Complete')

def readScore(scoreFile, ut, lt):
	scorelist = []
	with open(scoreFile) as f:
		print('Reading score file')
		for line in f:
			scorelist.append(float(line))
	up = np.percentile(scorelist, ut)  
	print('Top ', ut, ' percentile is: ', up)
	lp = np.percentile(scorelist, lt)  
	print('Bottom ', lt, ' percentile is: ', lp)
	print('Score Parse Complete')
	
def parseScoreNew(mirandaOut, outputFile, so, ut, lt):
	scorelist = []
	if so:
		with open(mirandaOut) as f, open(outputFile, "w") as text_file:
			print('Reading scores from miranda output file')
			for line in f:
				if line[0]=='>':
					lineSplit = line.split('\t')
					score = lineSplit[2]
					print('{}'.format(score), file=text_file)
		print("Success")
	else:
		with open(mirandaOut) as f:
			print('Processing miranda output file')
			for line in f:
				if line[0]=='>':
					lineSplit = line.split('\t')
					scorelist.append(float(lineSplit[2]))
		if outputFile!="off":
			try:
				with open(outputFile, "a") as text_file:
					print('Writing scores to output file')
					for score in scorelist:
						print('{}'.format(score), file=text_file)
			except:
				print('Could not print to output file')
		up = np.percentile(scorelist, ut)  
		print('Top ', ut, ' percentile is: ', up)
		lp = np.percentile(scorelist, lt)  
		print('Bottom ', lt, ' percentile is: ', lp)
		print('Score Parse Complete')