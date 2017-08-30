from setuptools import setup


setup(
	name='miRanda_input_processing',
	version='1.0',
	py_modules=['main'],
	install_requires=[
		'Click',
	],
	entry_points='''
		[console_scripts]
		main=main:cli
		hello=hello:cli
	''',
)
