# -*- Mode: Makefile -*- 

# the location of the Chombo "lib" directory
ifndef CHOMBO_HOME
    $(error Please define CHOMBO_HOME - see installation instructions.)
endif

# User needs to specify GRTRESNA_SOURCE location, e.g. using
# export GRTRESNA_SOURCE = ~/GRTresna/Source

# trace the chain of included makefiles
makefiles += releasedExamples_AMRPoisson_execVariableCoefficient

## Define the variables needed by Make.example

# Choose a method USE_CTTK or USE_CTTKHybrid
#cxxcppflags = -DUSE_CTTK
cxxcppflags = -DUSE_CTTKHybrid

# the base name(s) of the application(s) in this directory
ebase = Main_AxiDilatonBH

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTools BoxTools

# input file for 'run' target
INPUT = params.txt

# application-specific targets
src_dirs := $(GRTRESNA_SOURCE) \
            $(GRTRESNA_SOURCE)/simd \
            $(GRTRESNA_SOURCE)/utils \
			$(GRTRESNA_SOURCE)/Core \
            ../../Source/Core \
			$(GRTRESNA_SOURCE)/Matter \
            ../../Source/Matter \
            $(GRTRESNA_SOURCE)/Methods \
            $(GRTRESNA_SOURCE)/Tools \
            $(GRTRESNA_SOURCE)/Variables \
            $(GRTRESNA_SOURCE)/TaggingCriteria \
	        $(GRTRESNA_SOURCE)/Operator \
            $(GRTRESNA_SOURCE)/Operator/SolverOperator 

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example
