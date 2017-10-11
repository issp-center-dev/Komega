ACLOCAL_AMFLAGS         = -I config/m4
MAINTAINERCLEANFILES    = Makefile.in aclocal.m4 configure

EXTRA_DIST = doc
SUBDIRS = src test app/src

docs:
	cd doc/library/ja/;make latexpdfja; make html
	cd doc/library/en/;make latexpdfja; make html
	cd doc/software/ja/;make latexpdfja; make html
	cd doc/software/en/;make latexpdfja; make html

pngdocs:
	cd doc/library/ja/;sed -i -e "s/mathjax/pngmath/g" conf.py; make latexpdfja; make html
	cd doc/library/en/;sed -i -e "s/mathjax/pngmath/g" conf.py; make latexpdfja; make html
	cd doc/software/ja/;sed -i -e "s/mathjax/pngmath/g" conf.py; make latexpdfja; make html
	cd doc/software/en/;sed -i -e "s/mathjax/pngmath/g" conf.py; make latexpdfja; make html

docs-clean:
	cd doc/library/ja/; make clean; sed -i -e "s/pngmath/mathjax/g" conf.py
	cd doc/library/en/; make clean; sed -i -e "s/pngmath/mathjax/g" conf.py
	cd doc/software/ja/; make clean; sed -i -e "s/pngmath/mathjax/g" conf.py
	cd doc/software/en/; make clean; sed -i -e "s/pngmath/mathjax/g" conf.py
