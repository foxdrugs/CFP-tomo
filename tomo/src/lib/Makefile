DIRS =  OperatorLib \
		RaytraceLib \
		WavefrontLib \
		SIRTLib \


SHELL = /bin/sh

OUTPUT_LIBDIR=../../lib
$(shell mkdir -p $(OUTPUT_LIBDIR) >/dev/null)

all: compile
compile:
	@echo --------------------------------
	@echo Compiling libraries...
	@echo --------------------------------
	@for lib in $(DIRS); do cd $$lib; echo compiling $$lib; make -s; cd ..; done
clean:
	@echo --------------------------------
	@echo Cleaning libraries...
	@echo --------------------------------
	@for lib in $(DIRS); do cd $$lib; echo cleaning $$lib;  make clean; cd ..; done
