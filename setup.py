from setuptools import setup

setup(
	name='miRNA_SNP_processing_suite',
	version='1.0',
	py_modules=[
		'process_chromosome', 
		'parse_score', 
		'parse_mir', 
		'compressmir', 
		'procfirstrepass', 
		'first_pass'
	],
	install_requires=[
		'Click',
		'numpy'
	],
	entry_points='''
		[console_scripts]
		procchrom=process_chromosome:cli
		parsescore=parse_score:cli
		parsemir=parse_mir:cli
		compressmir=compressmir:cli
		procfirstrepass=procfirstrepass:cli
		firstpass=first_pass:cli
	''',
)
