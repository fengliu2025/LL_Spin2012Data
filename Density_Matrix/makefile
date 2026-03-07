CXX = g++
CXXFLAGS =  -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC



ROOT_INC = -I$(shell root-config --incdir)
ROOT_LIBS = $(shell root-config --libs)

CXXFLAGS += $(ROOT_INC) $(PYTHIA8_INC)
LIBS = $(ROOT_LIBS)

TARGET = MixEventAnalyzer

all: $(TARGET)


$(TARGET): main_Analysis.o
	$(CXX) -o $(TARGET) main_Analysis.o $(LIBS)

main_Analysis.o: main_Analysis.C
	$(CXX) $(CXXFLAGS) -c main_Analysis.C -o main_Analysis.o

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean


