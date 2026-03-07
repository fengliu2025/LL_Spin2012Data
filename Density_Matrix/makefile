CXX = g++
CXXFLAGS =  -stdlib=libc++ -pthread -std=c++17 -m64 -I/opt/homebrew/Cellar/root/6.34.08_1/include/root



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


