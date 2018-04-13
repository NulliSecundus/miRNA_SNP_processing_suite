import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console, such as a count of the number of scores processed''')
def cli(mirandafile, output, verbose):
	try:
		parsemir(mirandafile, output, verbose)
	except:
		print('Error parsing miranda output file')

def parsemir(mirandaFile, outputFile, v):
	print('Started Compressing Miranda File')
	
	output_container = []
	container = [] 
	
	with open(mirandaFile) as f:
		if v:
			print('Opened miranda file')
		
		count = 0
		for line in f:
			if line[0:2]=='>h':
				container.append(line.rstrip())
				count += 1
				if v: 
					if count%10000==0:
						print("Processed " + str(count) + " scores")
			elif line[0:2]=='>>':
				# end case
				# copy top score from from container
				if len(container)==1 :
					output_container.append(container[0])
				else:
					# determine top score line
					# print top score line to output file
					topLine = ""
					topScore = -1
					
					for line in container:
						splitline = line.split("\t")
						compare = float(splitline[2])
						if compare > topScore :
							topScore = compare
							topLine = line
					
					output_container.append(topLine)
				# reset the container
				container = [] 
					
	with open(outputFile, 'a') as out:
		if v:
			print('Writing to output file')
			
		for line in output_container:
			print('{}'.format(line), file=out)
