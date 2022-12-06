INCDIR = $(GARFIELD_INSTALL)/include

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -lGarfield `root-config --glibs` -lGeom -lgfortran -lm

all: StarTubeByAnsys StarTubeByneBEM

StarTubeByAnsys: StarTubeByAnsys.C 
	$(CXX) $(CFLAGS) StarTubeByAnsys.C
	$(CXX) `root-config --cflags` -o StarTubeByAnsys StarTubeByAnsys.o $(LDFLAGS)
	rm StarTubeByAnsys.o

StarTubeByneBEM: StarTubeByneBEM.C 
	$(CXX) $(CFLAGS) StarTubeByneBEM.C
	$(CXX) `root-config --cflags` -o StarTubeByneBEM StarTubeByneBEM.o $(LDFLAGS)
	rm StarTubeByneBEM.o

.PHONY: clean
clean:
	rm StarTubeByneBEM
	rm StarTubeByAnsys
