from setuptools import setup


setup(
	name='miRanda_input_processing',
	version='1.0',
	py_modules=['process_chromosome', 'parse_score', 'parse_mir'],
	install_requires=[
		'Click',
		'numpy'
	],
	entry_points='''
		[console_scripts]
		procchrom=process_chromosome:cli
		parsescore=parse_score:cli
		parsemir=parse_mir:cli
	''',
)
