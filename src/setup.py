# from setuptools import setup, Distribution
from setuptools import Distribution
from numpy.distutils.core import Extension, setup
import platform
import glob
import os


def readme():
    with open('../README.rst') as f:
        return f.read()


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False


# wrap OS-specific `SUEWS_driver` libs
sysname = platform.system()
if sysname == 'Windows':
    lib_name = 'SUEWS_driver.pyd'
elif sysname == 'Darwin':
    lib_name = 'SUEWS_driver.so'
elif sysname == 'Linux':
    lib_name = 'SUEWS_driver.so'

# load SUEWS Fortran source files
dir_f95 = 'suews/SUEWS-SourceCode'
target_f95 = [
    os.path.join(dir_f95, f)
    for f in
    ['suews_ctrl_const.f95',
     'suews_ctrl_driver.f95']]
all_f95 = glob.glob(os.path.join(dir_f95, '*.f95'))
exclude_f95 = [
    os.path.join(dir_f95, f)
    for f in
    ['suews_c_wrapper.f95',
     'suews_ctrl_sumin.f95',
     'suews_program.f95']
]
other_f95 = list(
    set(all_f95)
    - set(target_f95)
    - set(exclude_f95)
)
other_obj = [f.replace('.f95', '.o') for f in other_f95]
src_f95 = target_f95 + other_f95
for f in target_f95 + other_obj:
    print(f)

ext_modules = [
    Extension('supy.SUEWS_driver',
              target_f95,
              extra_f90_compile_args=['-O3 -cpp'],
              f2py_options=[
                  '--quiet',
                  ('-DF2PY_REPORT_ATEXIT' if sysname == 'Linux' else '')],
              extra_objects=other_obj,
              extra_link_args=[(''if sysname == 'Linux' else '-static')])]

setup(name='supy',
      version='0.3.1',
      description='the SUEWS model that speaks python',
      long_description=readme(),
      url='https://github.com/sunt05/SuPy',
      author='Ting Sun',
      author_email='ting.sun@reading.ac.uk',
      license='GPL-V3.0',
      packages=['supy'],
      package_data={
          'supy': [
              # lib_name,
           '*.json'
          ]
      },
      distclass=BinaryDistribution,
      ext_modules=ext_modules,
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'f90nml'
      ],
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
