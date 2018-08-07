import click
import numpy as np

@click.command()
@click.argument('mirandafile')
@click.option('-out', default='off', help='Print the parsed scores to the specified output file')
@click.option('--so', is_flag=True, help='Scores only option (for large files). Will only output scores without stats')
@click.option('--verbose', is_flag=True, help='Keep a count of the number of scores processed and print to console')
@click.option('-ut', default=90.0, help='The upper threshold in which the top percentile is calculated\nDefault 90.0')
@click.option('-lt', default=10.0, help='The lower threshold in which the bottom percentile is calculated\nDefault 10.0')
@click.option('--scorefile', is_flag=True, help='Flag if a pre-processed score file is provided instead of miranda output')
@click.option('--parsefile', is_flag=True, help='Flag if a parsed miranda output file is provided as input')
@click.option('-split', default=300000000, help='Splits scores output by the specified number per file\nDefault 300,000,000')
def cli(mirandafile, out, so, verbose, ut, lt, scorefile, parsefile, split):
	try:
		if scorefile:
			readScore(mirandafile, verbose, ut, lt)
		elif parsefile:
			parseScoreNew(mirandafile, out, so, verbose, ut, lt, split)
		else:
			parseScore(mirandafile, out, verbose, ut, lt)
		print('Success')
	except:
		print('Error')

def parseScore(mirandaOut, outputFile, v, ut, lt):
	scorelist = []
	with open(mirandaOut) as f:
		print('Processing miranda output file')
		count = 0
		for line in f:
			if line[0:2]=='>>':
				lineSplit = line.split('\t')
				scorelist.append(float(lineSplit[4]))
				if v:
					count += 1
					if count%10000==0:
						print(count)
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

def readScore(scoreFile, v, ut, lt):
	scorelist = []
	with open(scoreFile) as f:
		print('Reading score file')
		count = 0
		for line in f:
			scorelist.append(float(line))
			if v:
				count += 1
				if count%10000==0:
					print(count)
	up = np.percentile(scorelist, ut)  
	print('Top ', ut, ' percentile is: ', up)
	lp = np.percentile(scorelist, lt)  
	print('Bottom ', lt, ' percentile is: ', lp)
	print('Score Parse Complete')
	
def parseScoreNew(mirandaOut, outputFile, so, v, ut, lt):
	scorelist = []
	if so:
		count=0
		with open(mirandaOut) as f, open(outputFile, "w") as text_file:
			print('Reading scores from miranda output file')
			for line in f:
				if line[0]=='>':
					count+=1
					lineSplit = line.split('\t')
					score = lineSplit[2]
					print('{}'.format(score), file=text_file)
		print("Success")
	else:
		with open(mirandaOut) as f:
			print('Processing miranda output file')
			count = 0
			for line in f:
				if line[0]=='>':
					lineSplit = line.split('\t')
					scorelist.append(float(lineSplit[2]))
					if v:
						count += 1
						if count%10000==0:
							print(count)
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