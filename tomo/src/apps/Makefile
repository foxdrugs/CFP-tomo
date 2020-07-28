DIRS=Createoperators

SHELL = /bin/sh
all: compile
compile:
	@echo --------------------------------
	@echo Compiling apps...
	@echo --------------------------------
	@for lib in $(DIRS); do cd $$lib; make -s; make install; cd ..; done
clean:
	@echo --------------------------------
	@echo Cleaning apps...
	@echo --------------------------------
	@for lib in $(DIRS); do cd $$lib; make clean; cd ..; done
