from setuptools import setup, find_packages

with open('README.md') as readme_file:
    README = readme_file.read()

setup_args = dict(
    name='annotationsummary',
    version='0.0.1',
    description='Compute Summary and Comparison of MG-Rast Annotations and Diamond-Megan-Pipeline Annotations',
    long_description_content_type="text/markdown",
    long_description=README,
    license='MIT',
    packages=find_packages(),
    url='https://github.com/mariorauh/MasterThesis-PGPT/tree/master/Code/AS/annotationsummary',
    author='Mario Rauh',
    author_email='mario.rauh@student.uni-tuebingen.de',
    keywords=['KEGG','COG','Megan6','DIAMOND','Bioinformatics','AnnotationSummary', 'Genome Annotation'],
    zip_safe=False
)

install_requires = [
    'setuptools~=49.6.0',
    'pandas~=1.2.4',
    'matplotlib~=3.4.1',
    'numpy~=1.19.2',
    'scipy~=1.6.2'
]

if __name__ == '__main__':
    setup(**setup_args, install_requires=install_requires)
