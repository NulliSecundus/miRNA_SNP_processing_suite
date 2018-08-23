miRNA SNP Processing Suite

Copyright (C) (2018) Yang Laboratory, University of Wisconsin - Madison

See the file 'INSTALLATION' for installation details

Author: Brendan Drackley

Email: brendandrackley@gmail.com
 
Please send bug reports to: jyang75@wisc.edu

Primary Commands:

    procchrom [OPTIONS] REPORTFILE SNPFILE OUTPUT
		Arguments: 
			REPORTFILE (chr_N.txt) - INPUT chromosome report file containing SNP 
				reference information. 
				See ftp://ftp.ncbi.nlm.nih.gov/snp/00readme.txt 
			SNPFILE (rs_chrN.fas) - INPUT SNP sequence files by chromosome in fasta 
				format. See ftp://ftp.ncbi.nlm.nih.gov/snp/00readme.txt 
			OUTPUT (proc_chrN.fasta) - processed fasta formatted file 
				containing entries with sequences for each validated SNP. 
		
		Extracts validated SNPs with a unique genomic position located within 
		the 5' or 3' UTR or CDS of any gene. Trims rs_chrN.fas entries to 25 bp
		each side of SNP and creates separate FASTA entries for each allele.

    firstpass [OPTIONS] SNPFILE MIRNAFILE OUTPUT SCORE
		Arguments:
			SNPFILE - INPUT fasta-formatted processed chromosomal SNP file 
			MIRNAFILE - INPUT fasta-formatted mature miRNA file 
			OUTPUT - final processed file from first pass miRanda
			SCORE - score threshold parameter for miRanda 
			
		Wrapper interface for running a miranda first pass. For the given 
		snpFile, mirnaFile, and score threshold, iteratively runs miranda in 
		parallel on small subsections of the input files and combines the 
		outputs to a final processed output file. 
		
	procfirstpass [OPTIONS] MIRANDAFILE OUTPUT
		Arguments:
			MIRANDAFILE - INPUT processed final output from firstpass script
			OUTPUT - fully processed first pass miRanda output  
		
		Parses a miranda output file for a given chromosome and only retains 
		entries that have at least one, but not all, SNP alleles for a given 
		SNP ID in the upper 80th percentile. 

    secondpass [OPTIONS] MIRANDAFILE PROCSNPFILE MIRNAFILE OUTPUT
		Arguments:
			MIRANDAFILE - INPUT processed output from procfirstpass script 
			PROCSNPFILE - INPUT processed SNP file from procchrom script 
			MIRNAFILE - INPUT fasta-formatted mature miRNA file 
			OUTPUT - processed output from second pass of miRanda 
		
		Wrapper interface for running a miranda second pass.
		
	procsecondpass [OPTIONS] MIRANDAFILE OUTPUT
		Arguments:
			MIRANDAFILE - INPUT output file from secondpass script 
			OUTPUT - fully processed second pass miRanda output 
		
		Parses the miranda output for the second pass.

Utility Commands: 

    parsescore [OPTIONS] MIRANDAFILE
		Arguments:
			MIRANDAFILE - miRanda output file or other acceptable format file 
				containing miRanda score information. 
				
				Acceptable file formats are *
				* direct miRanda output containing '>>' summary lines
				* parsed miRanda output containing only data lines 
				* pre-processed score file containing only parsed scores 
		
		Utility command for generating statistical data from miranda output. 
		
		Unless specified in the options, the script will assume standard 
		miRanda output has been provided as the input file. 

		Options:
			-ut FLOAT	The upper threshold in which the top percentile is calculated.
						Default 80.0
			-lt FLOAT	The lower threshold in which the bottom percentile is calculated.
						Default 20.0
			-out TEXT	Print the parsed scores to the specified output file.
			--so		Scores only option (for large files). Will only output scores without stats.
			--scorefile	Flag if a pre-processed score file is provided instead of miranda output.
			--parsefile	Flag if a parsed miranda output file is provided as input. 
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
	
	
