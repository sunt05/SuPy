from setuptools import setup
from supy.version import __version__
# from setuptools import Distribution
# from numpy.distutils.core import setup
# import platform
# import glob
# import os


def readme():
    with open('../README.md', encoding='utf-8') as f:
        return f.read()


setup(name='supy',
      version=__version__,
      description='the SUEWS model that speaks python',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/sunt05/SuPy',
      author='Ting Sun',
      author_email='ting.sun@reading.ac.uk',
      license='GPL-V3.0',
      packages=['supy'],
      package_data={
          'supy':
          [
              'sample_run/*',
              'sample_run/Input/*',
              '*.json'
          ]
      },
      # distclass=BinaryDistribution,
      ext_modules=[],
      install_requires=[
          'numpy>=1.15.2',
          'pandas>=0.23.4',
          'scipy',
          'f90nml',
          'matplotlib',
          'seaborn',
          'supy_driver>=2018b20'  # a separate f2py-based driver
      ],
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      python_requires='~=3.5',
      classifiers=[
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
      ],
      zip_safe=False)
