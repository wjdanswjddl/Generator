# Make sure to set up GENIE before building
ROOT_CXX_FLAGS := $(shell root-config --cflags) -fPIC
ROOT_LIBS := $(shell root-config --libs --evelibs) -lMathMore

GENIE_CXX_FLAGS := -I$(shell genie-config --topsrcdir) -I${LOG4CPP_INC}
GENIE_LIBS := $(shell genie-config --libs) -L${LOG4CPP_LIB} -lxml2 -llog4cpp
GENIE_LIBS += -lEGPythia6 $(gsl-config --libs) -lGenVector
GENIE_LIBS += -L${LHAPDF_LIB} -lLHAPDF

# When including the Fortran code, you will likely need to link the executable
# to -lgfortran as well
test_xsec: test_xsec.o
	g++ $(ROOT_LIBS) $(GENIE_LIBS) -o $@ $^

%.o: %.cxx
	g++ $(ROOT_CXX_FLAGS) $(GENIE_CXX_FLAGS) -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.o test_xsec
