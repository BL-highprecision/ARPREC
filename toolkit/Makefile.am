if HAVE_FORTRAN
AM_CPPFLAGS=-I$(top_builddir)
AM_FFLAGS=@ax_f90_modflag@$(top_builddir)/fortran
LDADD = $(top_builddir)/fortran/libarprec_f_main.a \
        $(top_builddir)/fortran/libarprecmod.a \
        $(top_builddir)/src/libarprec.a $(FCLIBS)

TOOLKIT = mathinit$(EXEEXT) mathtool$(EXEEXT)
EXTRA_PROGRAMS = $(TOOLKIT)
CLEANFILES = $(TOOLKIT)

mathinit_SOURCES = globdata.f mathinit.f
nodist_mathinit_SOURCES = second.f
nodist_EXTRA_mathinit_SOURCES = dummy.cpp
mathinit_LINK = $(CXXLINK)

mathtool_SOURCES = mathtool.f globdata.f pslqsub.f quadsub.f \
                   rootsub.f zetapz.f 
nodist_mathtool_SOURCES = second.f
nodist_EXTRA_mathtool_SOURCES = dummy.cpp
mathtool_LINK = $(CXXLINK)

toolkit: $(TOOLKIT)

endif
