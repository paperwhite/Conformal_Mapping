MAKEFILE= Makefile

SRCS = ../Core/Edge.cpp ../Core/Face.cpp ../Core/FException.cpp ../Core/HalfEdge.cpp ../Core/OBJFileReader.cpp ../Core/Parser.cpp ../Core/Point.cpp ../Core/Solid.cpp ../Core/SolidDelegate.cpp ../Core/StringTokenizer.cpp ../Core/TopologyException.cpp ../Core/Trait.cpp ../Core/Vertex.cpp

INCS_PATHS = -I../Core

DRIVER = Tuette.cpp  
CMD = Tuette 

CC	= /usr/bin/g++ 

all:
	$(CC) -o $(CMD) $(INCS_PATHS) $(DRIVER) $(SRCS)

# DO NOT DELETE
