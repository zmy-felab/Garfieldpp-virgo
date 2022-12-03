INCDIR = $(GARFIELD_INSTALL)/include

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -lGarfield `root-config --glibs` -lGeom -lgfortran -lm

all: generate read
.PHONY: all
generate: generate.C 
	$(CXX) $(CFLAGS) generate.C
	$(CXX) `root-config --cflags` -o generate generate.o $(LDFLAGS)
	rm generate.o
read: read.C 
	$(CXX) $(CFLAGS) read.C
	$(CXX) `root-config --cflags` -o read read.o $(LDFLAGS)
	rm read.o
