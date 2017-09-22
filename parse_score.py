import click
import numpy as np

@click.command()
@click.argument('mirandafile')
@click.option('-o', default='off', help='Print the parsed scores to the specified output file')
@click.option('--verbose', is_flag=True, help='Keep a count of the number of scores processed and print to console')
@click.option('-ut', default=90.0, help='The upper threshold in which the top percentile is calculated')
@click.option('-lt', default=10.0, help='The lower threshold in which the bottom percentile is calculated')
def cli(mirandafile, o, verbose, ut, lt):
	try:
		parseScore(mirandafile, o, verbose, ut, lt)
	except:
		print('Error parsing scores')

def parseScore(mirandaOut, outputFile, v, ut, lt):
	scorelist = []
	with open(mirandaOut) as f:
		print('miranda output file opened')
		count = 0
		for line in f:
			if line[0]=='>':
				lineSplit = line.split('\t')
				scorelist.append(float(lineSplit[2]))
				if v:
					count += 1
					if count%1000==0:
						print(count)
	if outputFile!="off":
		try:
			with open(outputFile, "a") as text_file:
				print('output file opened')
				for score in scorelist:
					print('{}'.format(score), file=text_file)
		except:
			print('Could not print to output file')
    up = np.percentile(scorelist, ut)  
	print('Top ', ut, ' percentile is: ', up)
	lp = np.percentile(scorelist, lt)  
	print('Bottom ', lt, ' percentile is: ', lp)
	print('Score Parse Complete')