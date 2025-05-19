UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CXX = g++-14
  CXXFLAGS = -O3 -Wall -Wextra -std=c++17 -I ./toofus
  GLLINK = `pkg-config --libs gl glu glut`
	GLFLAGS = `pkg-config --cflags gl glu glut`	
	# on apple, use brew to install freeglut and mesa-glu
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++17 -I ./toofus
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif



# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif



# The compiler to be used
#CXX = g++-14

# The list of flags passed to the compiler
#CXXFLAGS = -O3 -Wall -Wextra -std=c++17 -I ~/toofus

# The list of source files
SOURCES = FiberCell2DSimulation.cpp Loading.cpp PeriodicCell.cpp \
Platelet.cpp InteractionPlatelet.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)

# All source files (listed in SOURCES) will be compiled into an object file 
# with the following command
#%.o:%.cpp
#	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean

all: FibresGrains

clean:
	@echo "\033[0;32m-> Remove object files\033[0m"
	rm -f *.o
	@echo "\033[0;32m-> Remove compiled applications\033[0m"
	rm -f FibresGrains

%.o: %.cpp 
	@echo "\033[0;32m-> COMPILING OBJECT" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
libFibresGrains.a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING LIBRARY" $@ "\033[0m"
	ar rcs $@ $^

FibresGrains: FibresGrains.cpp libFibresGrains.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o run.o
	$(CXX) -o $@ run.o libFibresGrains.a

see: see.cpp libFibresGrains.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFLAGS)
	$(CXX) -o $@ see.o libFibresGrains.a $(GLLINK)

build: build.cpp
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) $< -o $@
	