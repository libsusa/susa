CC = gcc
CXX = g++
CFLAG = -x c++ -c
CFLAGXX = -c -Wall
LIBFILES = rng.cpp base.cpp mt.cpp rrcosine.cpp convolutional.cpp utility.cpp modulation.cpp svd.cpp matrix.cpp
LIBOBJFILES = obj/rng.o obj/base.o obj/mt.o obj/rrcosine.o obj/convolutional.o obj/utility.o obj/modulation.o obj/svd.o obj/matrix.o

all :
	@echo . compile
	
	$(CXX) $(CFLAGXX) src/rng.cpp -o obj/rng.o
	$(CXX) $(CFLAGXX) src/base.cpp -o obj/base.o
	$(CXX) $(CFLAGXX) src/mt.cpp -o obj/mt.o
	$(CXX) $(CFLAGXX) src/rrcosine.cpp -o obj/rrcosine.o
	$(CXX) $(CFLAGXX) src/convolutional.cpp -o obj/convolutional.o
	$(CXX) $(CFLAGXX) src/utility.cpp -o obj/utility.o
	$(CXX) $(CFLAGXX) src/modulation.cpp -o obj/modulation.o
	$(CXX) $(CFLAGXX) src/svd.cpp -o obj/svd.o
	$(CXX) $(CFLAGXX) src/matrix.cpp -o obj/matrix.o


	@echo . archive

	ar rcs lib/libsusa.a $(LIBOBJFILES)

clean :
	rm obj/*.o
	rm lib/*.a
