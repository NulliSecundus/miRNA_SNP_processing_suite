import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
def cli(mirandafile, output):

	"""
	\b 
	Arguments:
		MIRANDAFILE - INPUT output file from secondpass script 
		OUTPUT - fully processed second pass miRanda output 
	
	Parses the miranda output for the second pass."""

	try:
		parsemir(mirandafile, output)
	except:
		print('Error parsing miranda output file')

def parsemir(mirandaOut, outputFile):
	print('Started ParseMir')
	
	# TODO
