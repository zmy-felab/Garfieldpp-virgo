INCDIR = $(GARFIELD_INSTALL)/include

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -lGarfield `root-config --glibs` -lGeom -lgfortran -lm

TARGET = StarTube

MODE = Ansys

$(TARGET)By$(MODE): $(TARGET)By$(MODE).C 
	$(CXX) $(CFLAGS) $(TARGET)By$(MODE).C
	$(CXX) `root-config --cflags` -o $(TARGET)By$(MODE) $(TARGET)By$(MODE).o $(LDFLAGS)
	rm $(TARGET)By$(MODE).o

MODE = neBEM

$(TARGET)By$(MODE): $(TARGET)By$(MODE).C 
	$(CXX) $(CFLAGS) $(TARGET)By$(MODE).C
	$(CXX) `root-config --cflags` -o $(TARGET)By$(MODE) $(TARGET)By$(MODE).o $(LDFLAGS)
	rm $(TARGET)By$(MODE).o

.PHONY: clean
clean:
	rm $(TARGET)*
