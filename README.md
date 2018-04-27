miRNA SNP Processing Suite V1.0

Primary Commands:

    procchrom [OPTIONS] REPORTFILE SNPFILE OUTPUT
		Inputs: REPORTFILE (chr_N.txt), SNPFILE (rs_chrN.fas)
		
        Extracts validated SNPs with a unique genomic position located within 
		the 5' or 3' UTR or CDS of any gene. Trims rs_chrN.fas entries to 25 bp
		each side of SNP and creates separate FASTA entries for each allele.

    parsemir [OPTIONS] MIRANDAFILE OUTPUT
        Compresses a miranda output file by only retaining the top scoring hit 
		from each miRNA-SNP pair and outputting to a text file. 
		
	procfirstpass
		Process the (parsed) first pass miranda output

    secondpass
        Process the second pass miranda output

Utility Commands: 

    parsescore [OPTIONS] MIRANDAFILE
	    Utility command for generating statistical data from miranda output.
		
		Options:
		  -so TEXT		Print the parsed scores to the specified output file
		  --verbose		Keep a count of the number of scores processed and 
						print to console
		  -ut FLOAT		The upper threshold in which the top percentile is 
						calculated
		  -lt FLOAT		The lower threshold in which the bottom percentile is 
						calculated
		  --help		Show this message and exit.

Files:

	compressmir.py - python source file for parsemir command
	
	parse_mir.py - [DEPRECIATED] python source file for parsemir command
	
	parse_score.py - python source file for parsescore command
	
	process_chromosome.py - python source file for procchrom command
	
	procmirout.py - python source file for procmirout command
	
	procsecondpass.py - python source file for secondpass command
	
	setup.py - setuptools python source file
	
	
