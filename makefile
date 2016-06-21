
all:lib test

lib:
	cd src;make

tests:
	cd test;make

clean:
	cd src;make clean
	cd test;make clean
