import click

@click.command()
@click.argument('mirandafile')
@click.argument('output')
@click.option('--verbose', is_flag=True, help='''Output additional information to
	the console, such as a count of the number of scores processed''')
def cli(mirandafile, output, verbose):
	# command line interface setup
	try:
		parsemir(mirandafile, output, verbose)
	except:
		print('Error parsing miranda output file')

def parsemir(mirandaFile, outputFile, v):
	print('Compressing Miranda File')
	
	output_container = [] # container for lines of final output
	container = [] # temporary container for comparison within a group
	
	with open(mirandaFile) as f:
		if v:
			print('Opened miranda file')
		
		count = 0
		sumcount = 0
		
		for line in f:
			# for each line in the miranda output file
			if line[0:2]=='>h':
				# if line is a data line
				container.append(line.rstrip())
				count += 1
				if v: 
					if count%10000==0:
						print("Processed " + str(count) + " scores")
						
			elif line[0:2]=='>>':
				# if line is a summary line
				# end case for a grouping
				# copy top score from from container
				sumcount += 1
				topLine = ""
				
				if len(container)==1 :
					# if the group only has one data line 
					topLine = str(container[0])
					#output_container.append(container[0])
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
	
	with open(outputFile, 'a') as out:
		if v:
			print('Writing to output file')
			
		for line in output_container:
			# write each topLine to the output file 
			print('{}'.format(line), file=out)
