import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'kamodo',
    version = '0.0.1',
    author = 'Asher Pembroke',
    author_email = 'apembroke@gmail.com',
    description = 'A functional api for scientific data',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = 'https://ccmc.gsfc.nasa.gov/Kamodo/',
    packages = setuptools.find_packages(),
    classifiers = [
		"Programming Language :: Python :: 3.7",
	    "Operating System :: OS Independent",
	    "License :: OSI Approved",
    ],
    install_requires = [
		'numpy',
		'scipy',
		'sympy',
		'pandas',
		'plotly',
		'pytest',
		'psutil',
		'antlr4-python3-runtime',
    	],
    license='NASA OPEN SOURCE AGREEMENT VERSION 1.3',

)

