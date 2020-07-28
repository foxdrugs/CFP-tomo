# Makefile

CPP = mpic++

CPPFLAGS=-std=c++0x -g -Wall -m64 -O3 -fopenmp 

INCLUDEPATH=-I. \
			-I../../lib/OperatorLib \
			-I../../lib/RaytraceLib \
			-I../../lib/SIRTLib  \
			-I../../lib/WavefrontLib  

LIBSPATH= 	-L../../../lib 

LIBS= -lOperator \
	  -lRaytrace  \
	  -lWavefront \
	  -lSIRT

BINDIR=../../../bin

TARGET=test

DEPS = $(shell find ./ -name "*.h")
SRCS = $(shell find ./ -name "*.cpp")
OBJS = $(SRCS:%.cpp=%.o)


.SUFFIXES : .o .cpp
.cpp.o	:
	$(CPP) -c $^ $(INCLUDEPATH) $(CPPFLAGS)

all: $(TARGET)
$(TARGET):$(OBJS) $(DEPS)
	$(info ************Making APP************************)
	$(CPP) -o $@ $(OBJS) $(INCLUDEPATH) $(LIBSPATH) $(LIBS)  $(CPPFLAGS)
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

