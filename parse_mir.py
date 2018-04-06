import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
def cli(mirandafile, output):
	try:
		parsemir(mirandafile, output)
	except:
		print('Error parsing miranda output file')

def parsemir(mirandaOut, outputFile):
	print('Started ParseMir')
	with open(mirandaOut) as f:
		print('opened miranda output file')
		with open(outputFile, 'a') as out:
			print('opened output file')
			count = 0
			for line in f:
				if line[0]=='>':
					print('{}'.format(line), file=out)
					count += 1
					if count%10000==0:
						print(count)
