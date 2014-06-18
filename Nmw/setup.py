from setuptools import setup, find_packages
setup(
        name='nmw',
        version='0.1',
        author='Janina Mass',
        author_email='janina.mass@hhu.de',
        install_requires=["pickle"],
        packages=find_packages(),
        scripts=['nmw/nmw.py'],
        license='GPLv3',
        description='Simple Needleman-Wunsch algorithm implementation',
        long_description=open('README.txt').read(),
        classifiers=[
            'Topic :: Scientific/Engineering :: Bio-Informatics'
            ],
        )
