# -*- makefile -*-
.PHONY: main clean test pip supy

# OS-specific configurations
ifeq ($(OS),Windows_NT)
	PYTHON_exe = /c/Users/sunt05/Anaconda2/python.exe
	# F2PY_PY= /c/Users/sunt05/Anaconda2/Scripts/f2py.py
	# F2PY_EXE = $(PYTHON) $(F2PY_PY)
	TARGET=$(MODULE).pyd
else
	UNAME_S := $(shell uname -s)
	TARGET=$(MODULE).so

	ifeq ($(UNAME_S),Linux) # Linux
		PYTHON_exe=python
		# F2PY_EXE = f2py
	endif

	ifeq ($(UNAME_S),Darwin) # macOS
		PYTHON_exe=python
		# F2PY_EXE = f2py
	endif

endif

MODULE=SUEWS_driver

SUEWS_dir = suews/SUEWS-SourceCode

PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)
# All the files which include modules used by other modules (these therefore
# need to be compiled first)
FILES = SUEWS_const.f95  \
				SUEWS_driver.f95


# make fortran exe
main:
	$(MAKE) -C $(SUEWS_dir) clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) main; # make SUEWS with the `main` recipe
	-rm -rf *.o *.mod *.f95 *.a *.dSYM

# make fortran exe and run test cases
check:
	$(MAKE) -C $(SUEWS_dir) clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) check; # make SUEWS with the `main` recipe
	-rm -rf *.o *.mod *.f95 *.a *.dSYM

# make supy dist
supy:
	$(info $$PYTHON is [${PYTHON}])
	$(MAKE) -C $(SUEWS_dir) clean; # clean Fortran SUEWS build
	$(MAKE) -C $(SUEWS_dir) main; # make SUEWS with the `main` recipe
	$(PYTHON) setup.py bdist_wheel # all f2py compilation is done by `setup.py`
	-rm -rf *.o *.mod *.f95 *.a *.dSYM

# If wanted, clean all *.o files after build
clean:
	$(MAKE) -C $(SUEWS_dir) clean;
	 -rm -rf *.o *.mod *.dSYM $(TARGET) supy/$(MODULE).*;
	 -$(PYTHON) setup.py clean --all

# clean all existing builds, rebuild f2py libs, build wheels and test
test-supy:
	$(MAKE) clean;
	$(MAKE) supy;
	$(PYTHON) setup.py test

# clean all existing builds, rebuild f2py libs, build wheels and submit
pip:
	$(MAKE) clean;
	$(MAKE) supy;
	$(PYTHON) setup.py bdist_wheel upload
