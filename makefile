
all:lib tests

lib:
	cd src;make

tests:
	cd test;make

clean:
	cd src;make clean
	cd test;make clean
