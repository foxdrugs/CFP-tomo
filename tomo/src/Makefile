DIRS =  lib apps 

SHELL = /bin/sh
all: compile
compile:
	@echo --------------------------------
	@echo Compiling src...
	@echo --------------------------------
	@for lib in $(DIRS); do cd $$lib; echo compiling $$lib; make -s; cd ..; done
clean:
	@echo --------------------------------
	@echo Cleaning src...
	@echo --------------------------------
	@for lib in $(DIRS); do cd $$lib; echo cleaning $$lib; make clean; cd ..; done
