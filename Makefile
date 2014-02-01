CC = gcc
CXX = g++
CFLAG = -x c++ -c
ODIR = obj
SDIR = src
LDIR = lib
CFLAGXX = -c -Wall
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(addprefix $(ODIR)/,$(notdir $(SOURCES:%.cpp=%.o)))


$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(CFLAGXX) -o $@ $<

$(ODIR):
	mkdir -p $(ODIR)

$(LDIR):
	mkdir -p $(LDIR)

all: $(LDIR) $(ODIR) $(OBJECTS)
	ar rcs lib/libsusa.a $(OBJECTS)


clean:
	rm -f $(ODIR)/*.o
	rm -f $(LDIR)/*.a
	rm -d $(ODIR)
	rm -d $(LDIR)
