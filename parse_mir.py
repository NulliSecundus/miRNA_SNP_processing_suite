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
	
	container = [] 
	
	with open(mirandaOut) as f:
		print('opened miranda output file')
		with open(outputFile, 'a') as out:
			print('opened output file')
			count = 0
			for line in f:
				if line[0:2]=='>h':
					container.append(line)
					#print('{}'.format(line), file=out)
					count += 1
					if count%1000==0:
						print(count)
				elif line[0:2]=='>>':
					# end case
					# print top score from from container
					if len(container)==1 :
						print('{}'.format(container[0]), file=out)
					else:
						# determine top score line
						# print top score line to output file
						print("TODO") 
					# since end case reached
					# reset the container
					container = [] 
