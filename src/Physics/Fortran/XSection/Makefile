#
# Makefile for GENIE interface for Fortran-based cross section models
#
# Author: Steven Gardiner <gardiner \at fnal.gov>
#

# Make "all" the default target, even though we'll define "fortlib"
# before it
.PHONY: default
default: all

# Use gfortran as the default Fortran compiler
FC = gfortran

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

FORTRAN_SOURCES  := $(wildcard *.f90)
FORTRAN_OBJECTS  := $(addsuffix .$(ObjSuf), $(basename $(FORTRAN_SOURCES)))
FORTRAN_LIB_NAME := GFortranXS
FORTRAN_LIB_FILE := lib$(FORTRAN_LIB_NAME).$(DllSuf)
FORTRAN_LIB_FILE_WITH_PATH := $(GENIE_LIB_PATH)/$(FORTRAN_LIB_FILE)
#FORTRAN_LIB_FILE_WITH_IPATH := $(GENIE_LIB_INSTALLATION_PATH)/$(FORTRAN_LIB_FILE)

PACKAGE           = Physics/Fortran/XSection
PACKAGE_ABBREV    = PhFortranXS
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    = -l$(FORTRAN_LIB_NAME) -L$(GENIE_LIB_PATH)

fortobjs: $(FORTRAN_SOURCES)
	$(FC) -fPIC -c $(FORTRAN_SOURCES)

fortlib: fortobjs
	$(FC) -shared -fPIC -o $(FORTRAN_LIB_FILE_WITH_PATH) $(FORTRAN_OBJECTS)

all     : fortlib rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

FORCE:
