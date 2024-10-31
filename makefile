test : OcTree.o  test_FSI.o 
	g++ -o test OcTree.o test_FSI.o POINT2.o POINT3.o MarchingCubes.o Display.o

OcTree.o : OcTree.cpp
	g++ -O3 -c OcTree.cpp

test_FSI.o: test_FSI.cpp
	g++ -O3 -c test_FSI.cpp
