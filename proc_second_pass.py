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
	
	# TODO
