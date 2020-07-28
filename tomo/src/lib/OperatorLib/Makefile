# Makefile

CPP = g++
CPPFLAGS=-std=c++0x -Wall -m64 -O3 -fopenmp 
INCLUDEPATH=-I.  \
			-I../OperatorLib \
			-I../PerturbLib \
			-I../RaytraceLib \
			-I../SIRTLib \
			-I../WavefrontLib 

M_LIBS=../../../lib/libOperator.a

DEPS = $(shell find ./ -name "*.h")
SRCS = $(shell find ./ -name "*.cpp")

OBJS = $(SRCS:%.cpp=%.o)

.SUFFIXES : .o .cpp
.cpp.o	:
	$(CPP) -c $^ $(INCLUDEPATH) $(CPPFLAGS)

all: $(M_LIBS)
$(M_LIBS):$(OBJS) $(DEPS)
	$(info ************Making lib OperatorLib************************)
	@mkdir -p ../../../lib
	ar -rcs $(M_LIBS) $(OBJS)
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
	-rm -f $(M_LIBS)

