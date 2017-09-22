from setuptools import setup


setup(
	name='miRanda_input_processing',
	version='1.0',
	py_modules=['process_chromosome', 'parsescore'],
	install_requires=[
		'Click',
	],
	entry_points='''
		[console_scripts]
		procchrom=process_chromosome:cli
		parsescore=parse_score:cli
	''',
)
