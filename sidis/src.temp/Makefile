#-----------------------------------------------------------------------
# makefile for compiling sidis simulation code
#-----------------------------------------------------------------------
OBJECTS = main.o
CXXFLAGS = -std=c++17

denoise: $(OBJECTS)
	$(CXX) -o sidis.exe $(OBJECTS)

clean:
	@echo 'remove build files......'
	@rm -rf *.o *~ *.exe

%.o: %.cpp
	$(CXX) -fPIC $(CXXFLAGS) -c $< -O3
