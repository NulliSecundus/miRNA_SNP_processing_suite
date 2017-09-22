import click

@click.command()
@click.argument('mirandafile')
@click.argument('outputfile')
def cli(mirandafile, outputfile):
	try:
		parseScore(mirandafile, outputfile)
	except:
		print('Error parsing scores')

def parseScore(mirandaOut, outputFile):
	count = 0
	with open(mirandaOut) as f:
		print("miranda output file opened")
		with open(outputFile, "a") as text_file:
			print("output file opened")
			for line in f:
				if line[0]==">":
					lineSplit = line.split("\t")
					print("{}".format(lineSplit[2]), file=text_file)
					count += 1
					if count%1000==0:
						print(count)
	print("Score Parse Complete")

if __name__ == '__main__':
    pass

