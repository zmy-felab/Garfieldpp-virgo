INCDIR = $(GARFIELD_INSTALL)/include

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -lGarfield `root-config --glibs` -lGeom -lgfortran -lm

TARGET = CeramicGEM

$(TARGET): $(TARGET).C 
	$(CXX) $(CFLAGS) $(TARGET).C
	$(CXX) `root-config --cflags` -o $(TARGET) $(TARGET).o $(LDFLAGS)
	rm $(TARGET).o
.PHONY: clean
clean:
	rm $(TARGET)
