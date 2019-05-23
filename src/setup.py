from setuptools import setup
import pandas as pd

ser_ver = pd.read_json('./supy/supy_version.json', typ='series')
print(ser_ver)
__version__ = f'{ser_ver.ver_milestone}.{ser_ver.ver_major}.{ser_ver.ver_minor}{ser_ver.ver_remark}'


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
          'dask[complete]',  # needs all dask and its dependencies
          'f90nml',
          'matplotlib',
          'seaborn',
          'atmosp',  # my own `atmosp` module forked from `atmos-python`
          'cdsapi',
          'xarray',
          'click', # cmd tool
          'supy_driver>=2018rc8'  # a separate f2py-based driver
      ],
      entry_points={
          'console_scripts':[
              'suews-run=supy.cmd.SUEWS:SUEWS',
              'suews-convert=supy.cmd.table_converter:convert_table_cmd',
              ]},
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      python_requires='~=3.6',
      classifiers=[
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
      ],
      zip_safe=False)
