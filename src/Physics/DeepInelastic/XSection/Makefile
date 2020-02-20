#
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

PACKAGE           = Physics/DeepInelastic/XSection
PACKAGE_ABBREV    = PhDISXS
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

# Add extra source files (from Huma) to this library.
# We get them in this way to avoid problems generating ROOT
# dictionaries from her code.
EXTRA_SOURCES := $(wildcard *.cpp)
EXTRA_OBJECTS = $(addsuffix .$(ObjSuf), $(basename $(EXTRA_SOURCES)))

LIB_DEPEND := $(LIB_DEPEND) $(EXTRA_OBJECTS)

FORCE:
