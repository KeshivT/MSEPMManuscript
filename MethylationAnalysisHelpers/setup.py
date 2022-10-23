#! /usr/bin/env python3


try:
    from setuptools import setup
except ImportError as e:
    print('Please install python setuptools, '
          'https://packaging.python.org/tutorials/installing-packages/#use-pip-for-installing ')
    raise e


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='MethAnalysisHelpers',
      version='0.0.2',
      description='Analysis Processing Helpers',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/NuttyLogic/EPMTraitAssocition',
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='MIT',
      packages=['MethAnalysisHelpers',
                'MethAnalysisHelpers.NormalizationRefs'],
      classifiers=['Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8'],
      platforms=["Linux", "Mac OS-X", "Unix"],
      requires=['numpy', 'tqdm', 'joblib', 'scipy', 'seaborn', 'matplotlib'],
      install_requires=['numpy>=1.16.3', 'tqdm>=4.31.1',
                        'setuptools>=46.0.0', 'joblib>=0.14.0',
                        'seaborn>=0.11.0', 'scipy>=1.4.1',
                        'matplotlib>=3.3.2'],
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True,
      zip_safe=True
      )
