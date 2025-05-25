ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src tests include fortran toolkit

bin_SCRIPTS=arprec-config
docdir=${datadir}/doc/${PACKAGE}
dist_doc_DATA = README doc/arprec.pdf

cpp-demo:
	cd tests && make demo

if HAVE_FORTRAN

fortran-demo:
	cd fortran && make demo

demo: fortran-demo cpp-demo

toolkit:
	cd toolkit && make toolkit

else

fortran-demo:
	@echo "You need a Fortran 95 compiler to build fortran demo programs." && false

demo: cpp-demo

toolkit:
	@echo "You need a Fortran 95 compiler to build the toolkit." && false

endif

time:
	cd tests && make time

changelog:
	git log >ChangeLog

.PHONY: toolkit demo fortran-demo cpp-demo time changelog

