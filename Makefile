# -*- makefile -*-
.PHONY: main clean test pip supy docs

# OS-specific configurations
ifeq ($(OS),Windows_NT)
	PYTHON_exe = python.exe

else
	UNAME_S := $(shell uname -s)


	ifeq ($(UNAME_S),Linux) # Linux
		PYTHON_exe=python

	endif

	ifeq ($(UNAME_S),Darwin) # macOS
		PYTHON_exe=python

	endif

endif

src_dir = src
docs_dir = docs


PYTHON := $(if $(PYTHON_exe),$(PYTHON_exe),python)
# All the files which include modules used by other modules (these therefore
# need to be compiled first)

MODULE = supy

# default make options
main:
	$(MAKE) -C $(src_dir) main
	$(MAKE) -C $(docs_dir) html

# house cleaning
clean:
	$(MAKE) -C $(src_dir) clean
	$(MAKE) -C $(docs_dir) clean

# make supy and run test cases
supy:
	$(MAKE) -C $(src_dir) test

# make docs and open index
docs:
	$(MAKE) -C $(docs_dir) html
	open $(docs_dir)/build/html/index.html

# upload wheels to pypi using twine
upload:
	$(MAKE) -C $(src_dir) upload

# upload wheels to pypi using twine
livehtml:
	$(MAKE) -C $(docs_dir) livehtml