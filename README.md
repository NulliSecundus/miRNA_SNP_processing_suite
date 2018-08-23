miRNA SNP Processing Suite

Copyright (C) (2018) Yang Laboratory, University of Wisconsin - Madison

See the file 'INSTALLATION' for installation details

Author: Brendan Drackley

Email: brendandrackley@gmail.com
 
Please send bug reports to: jyang75@wisc.edu

Primary Commands:

    procchrom [OPTIONS] REPORTFILE SNPFILE OUTPUT
		Inputs: REPORTFILE (chr_N.txt), SNPFILE (rs_chrN.fas)
		
		Extracts validated SNPs with a unique genomic position located within 
		the 5' or 3' UTR or CDS of any gene. Trims rs_chrN.fas entries to 25 bp
		each side of SNP and creates separate FASTA entries for each allele.

    firstpass [OPTIONS] SNPFILE MIRNAFILE OUTPUT SCORE
		Wrapper interface for running a miranda first pass. For the given 
		snpFile, mirnaFile, and score threshold, iteratively runs miranda in 
		parallel on small subsections of the input files and combines the 
		outputs to a final processed output file. 
		
	procfirstpass [OPTIONS] MIRANDAFILE
		Parses a miranda output file for a given chromosome and only retains 
		entries that have at least one, but not all, SNP alleles for a given 
		SNP ID in the upper 80th percentile. 

    secondpass [OPTIONS] 
		Wrapper interface for running a miranda second pass.
		
	procsecondpass [OPTIONS]
		Parses the miranda output for the second pass. 

Utility Commands: 

    parsescore [OPTIONS] MIRANDAFILE
		Utility command for generating statistical data from miranda output.

		Options:
			-so TEXT	Print the parsed scores to the specified output file
			--verbose	Keep a count of the number of scores processed and 
						print to console
			-ut FLOAT	The upper threshold in which the top percentile is 
						calculated
			-lt FLOAT	The lower threshold in which the bottom percentile is 
						calculated
			--help		Show this message and exit.

Files:

	first_pass.py - python source file for firstpass command
	
	parse_score.py - python source file for parsescore command
	
	proc_first_pass.py - python source file for procfirstpass command
	
	proc_first_repass.py - python source file for procfirstrepass command
	
	proc_second_pass.py - python source file for procsecondpass command
	
	process_chromosome.py - python source file for procchrom command
	
	second_pass.py - python source file for secondpass command
	
	setup.py - setuptools python source file
	___________________________________________________________________________
	
	parse_mir.py - [DEPRECIATED] python source file for parsemir command
	
	compress_mir.py - [DEPRECIATED] python source file for compressmir command
	
	
