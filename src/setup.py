from setuptools import setup
# from setuptools import Distribution
# from numpy.distutils.core import setup
# import platform
# import glob
# import os


def readme():
    with open('../README.md') as f:
        return f.read()


setup(name='supy',
      version='0.5.2',
      description='the SUEWS model that speaks python',
      long_description=readme(),
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
          'numpy',
          'pandas',
          'scipy',
          'f90nml',
          'suews_driver>=2018b6'  # a separate f2py-based driver
      ],
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
