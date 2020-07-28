# Makefile

CPP = g++

CPPFLAGS=-std=c++0x -Wall -m64 -O3 -fopenmp 

INCLUDEPATH=-I. \
			-I../Base \
			-I/opt/OSF/EIGEN \

LIBSPATH=

LIBS=  

BINDIR=./

TARGET=Exe_Point

DEPS = $(shell find ./ -name "*.h")
SRCS = $(shell find ./ -name "*.cpp")
OBJS = $(SRCS:%.cpp=%.o)

.SUFFIXES : .o .cpp
.cpp.o	:
	$(CPP) -c $^ $(INCLUDEPATH) $(CPPFLAGS)

all: $(TARGET)
$(TARGET):$(OBJS) $(DEPS)
	$(info ************Making APP ************************)
	$(CPP) -o $@ $(OBJS) $(INCLUDEPATH) $(LIBS) $(LIBSPATH) $(CPPFLAGS)

	@echo "          "
	@echo "source files:"$(SRCS)
	@echo "          "
	@echo "object files:"$(OBJS)
	@echo "          "
	@echo "include files:"$(DEPS)
	@echo "          "

.PHONY:clean
clean:
	-rm -f $(OBJS)
	-rm -f $(BINDIR)/$(TARGET)

install:
	mv $(TARGET) $(BINDIR)/$(TARGET)


