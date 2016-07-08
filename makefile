
all:lib tests apps

lib:
	cd src/;make

apps:
	cd app/src/;make

tests:
	cd test/;make

clean:
	cd src/;make clean
	cd test/;make clean
	cd app/src/;make clean
