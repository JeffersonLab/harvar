#-----------------------------------------------------------------------
# makefile for compiling sidis simulation code
#-----------------------------------------------------------------------
OBJECTS = main.o
CXXFLAGS = -std=c++17

all:
	g++ -O2 -std=c++17 main.cpp -I. `root-config --cflags --glibs` -o generator.exe
clean:
	@rm -rf *.o generator.exe
