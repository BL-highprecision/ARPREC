LDADD = $(top_builddir)/src/libarprec.a
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/include -I$(top_srcdir)/include

TESTS = mp_test c_mp_test io
check_PROGRAMS = mp_test c_mp_test io trig
EXTRA_PROGRAMS = pslq1 pslq2 pslq3 mpslq1 mpslq2 mpslq3 quads mp_timer mp_ex
CLEANFILES = $(EXTRA_PROGRAMS) test

mp_test_SOURCES = mp_test.cpp

c_mp_test_SOURCES = c_mp_test.c
c_mp_test_LINK = $(CXXLINK)

if HAVE_FORTRAN
TESTS += f_mp_test
check_PROGRAMS += f_mp_test
AM_FCFLAGS = @ax_f90_modflag@$(top_builddir)/fortran
AM_FFLAGS = $(AM_FCFLAGS)

f_mp_test_SOURCES = f_mp_test.f
f_mp_test_LINK = $(CXXLINK)
f_mp_test_LDADD = $(top_builddir)/fortran/libarprecmod.a \
               $(top_builddir)/fortran/libarprec_f_main.a \
               $(LDADD) $(top_builddir)/src/libarprec.a $(FCLIBS)
endif

mp_timer_SOURCES = mp_timer.cpp tictoc.cpp

io_SOURCES = io.cpp
trig_SOURCES = trig.cpp
mp_ex_SOURCES = mp_ex.cpp
pslq1_SOURCES = pslq1.cpp pslq1_main.cpp pslq_main_util.cpp pslq1_util.cpp \
                tictoc.cpp
pslq2_SOURCES = pslq2.cpp pslq2_main.cpp pslq_main_util.cpp pslq1_util.cpp \
                pslq2_util.cpp tictoc.cpp
pslq3_SOURCES = pslq3.cpp pslq3_main.cpp pslq_main_util.cpp pslq1_util.cpp \
                pslq2_util.cpp tictoc.cpp
mpslq1_SOURCES = mpslq1.cpp mpslq1_main.cpp pslq_main_util.cpp \
                 mpslq1_util.cpp pslq1_util.cpp tictoc.cpp
mpslq2_SOURCES = mpslq2.cpp mpslq2_main.cpp pslq_main_util.cpp \
                 mpslq1_util.cpp pslq1_util.cpp pslq2_util.cpp tictoc.cpp
mpslq3_SOURCES = mpslq3.cpp mpslq3_main.cpp pslq_main_util.cpp \
                 mpslq1_util.cpp pslq1_util.cpp pslq2_util.cpp tictoc.cpp
quads_SOURCES = quads.cpp util.cpp arprec-integrate.cpp quad-erf.cpp \
                quad-gs.cpp quad-ts.cpp

noinst_HEADERS = arprec-integrate.h integrate.h matrix.h mpslq1.h mpslq2.h \
		 mpslq3.h pslq1.h pslq2.h pslq3.h pslq_main.h quad-erf.h \
		 quad-gs.h quad-ts.h tictoc.h util.h pslq1_templates.cpp \
		 pslq2_templates.cpp pslq3_templates.cpp mpslq1_templates.cpp \
		 mpslq2_templates.cpp mpslq3_templates.cpp

demo: $(EXTRA_PROGRAMS)

time: mp_timer$(EXEEXT)
	./mp_timer$(EXEEXT)

.PHONY: time demo

